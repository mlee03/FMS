
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
<<<<<<< HEAD
#include <fms_platform.h>
=======
>>>>>>> b402a7097b2ec57cf3b0aafff80ccfad4773a20f

  use mpp_mod, only : mpp_init, mpp_pe, mpp_npes, mpp_root_pe
  use mpp_mod, only : mpp_sync
  use mpp_mod, only : mpp_set_stack_size, mpp_init_test_requests_allocated
<<<<<<< HEAD
  use mpp_mod, only : mpp_sum, mpp_exit
  use mpp_mod, only : mpp_error, FATAL

  implicit none

  integer                                      :: n
  real(FLOAT_KIND), allocatable, dimension(:)  :: a4, b4, c4
  real(DOUBLE_KIND), allocatable, dimension(:) :: a8, b8, c8
  integer                                      :: i, pe, npes, root
  integer, allocatable, dimension(:)           :: pelist
=======
  use mpp_mod, only : mpp_sum
  use mpp_mod, only : mpp_error, FATAL
  use platform_mod

  implicit none

  integer                                       :: ierr
  integer                                       :: i, pe, npes, root, fullsum, pesum
  integer, allocatable, dimension(:)            :: pelist
>>>>>>> b402a7097b2ec57cf3b0aafff80ccfad4773a20f

  call mpp_init(mpp_init_test_requests_allocated)
  call mpp_set_stack_size(3145746)
  pe = mpp_pe()
  npes = mpp_npes()
  root = mpp_root_pe()
<<<<<<< HEAD

  if( pe.EQ.root ) print *, '------------------> Calling test_mpp_sum <------------------'
    call test_mpp_sum_dim1()
    call test_mpp_sum_pelist()
    call test_mpp_sum_large()
  if( pe.EQ.root ) print *, '------------------> Finished test_mpp_sum <------------------'

  call mpp_exit()

contains

  subroutine test_mpp_sum_dim1()

! Verify that mpp_sum calculates the sum for array with dim=1
  n=1
  allocate( a4(n), a8(n), b4(npes), b8(npes), c4(n), c8(n))

  a4 = real(pe+1, kind=FLOAT_KIND)
  a8 = real(pe+1, kind=DOUBLE_KIND)
  call mpp_sync()
  call mpp_sum(a4(1:n),n)
  call mpp_sum(a8(1:n),n)
  if (pe .EQ. root) then
    b4 = (/(i, i=1,npes, 1)/)
    b8 = (/(i, i=1,npes, 1)/)
    c4 = sum(b4)
    c8 = sum(b8)

    if (a4(1) .ne. c4(1)) call mpp_error(FATAL, "Real4: mpp_sum differs from fortran intrinsic sum")
    if (a8(1) .ne. c8(1)) call mpp_error(FATAL, "Real8: mpp_sum differs from fortran intrinsic sum")
  endif
  deallocate(a4, b4, c4, a8, b8, c8)

  end subroutine test_mpp_sum_dim1

  subroutine test_mpp_sum_pelist()

! Verify mpp_sum accurately sums from select list of pes
  allocate( a4(n), a8(n), b4(npes), b8(npes), c4(n), c8(n))
  allocate( pelist(0:npes-1) )

  a4 = real(pe+1, kind=FLOAT_KIND)
  a8 = real(pe+1, kind=DOUBLE_KIND)
  call mpp_sync()
  if (pe .LE. npes-2) then
    pelist = (/(i,i=0,npes-2)/)
    if (pe .EQ. root) print *,'PE: ',mpp_pe(),' pelist: ',pelist(:)
    call mpp_sum(a4(1:n),n, pelist=pelist)
    call mpp_sum(a8(1:n),n, pelist=pelist)
  endif
  if (pe .EQ. root) then
    b4 = (/(i, i=1,npes-1, 1)/)
    b8 = (/(i, i=1,npes-1, 1)/)
    c4 = sum(b4)
    c8 = sum(b8)
    print *, 'b4: ',b4, 'a4: ',a4, 'c4: ',c4
    if (a4(1) .ne. c4(1)) call mpp_error(FATAL, "Real4 with pelist: mpp_sum differs from fortran intrinsic sum")
    if (a8(1) .ne. c8(1)) call mpp_error(FATAL, "Real8 with pelist: mpp_sum differs from fortran intrinsic sum")
  endif
  deallocate(a4, b4, c4, a8, b8, c8)

  end subroutine test_mpp_sum_pelist

  subroutine test_mpp_sum_large()

! Verify mpp_sum works for large array dim with random values
  n=30000
  allocate( a4(n), a8(n) )

  call random_number(a4)
  call random_number(a8)

  call mpp_sync()
  call mpp_sum(a4(1:n),n)
  call mpp_sum(a8(1:n),n)

! Verify mpp_sum works for large array dim with random values and mpp_sum length less than array dim
  call random_number(a4)
  call random_number(a8)

  call mpp_sync()
  call mpp_sum(a4(1:n),n/2)
  call mpp_sum(a8(1:n),n/2)

  end subroutine test_mpp_sum_large
=======
  allocate( pelist(0:npes-2) )
  if (pe .LE. npes-2) pelist = (/(i,i=0,npes-2)/)
  if (pe .EQ. root) then
    fullsum = sum( (/(i, i=1,npes, 1)/) )
    pesum = sum( (/(i, i=1,npes-1, 1)/) )
  endif

  if( pe.EQ.root ) print *, '------------------> Calling test_mpp_sum <------------------'
    call test_mpp_sum_scalar(pe,npes,root,pelist,fullsum,pesum)
    call test_mpp_sum_2D(pe,npes,root,pelist,fullsum,pesum)
    call test_mpp_sum_3D(pe,npes,root,pelist,fullsum,pesum)
    call test_mpp_sum_4D(pe,npes,root,pelist,fullsum,pesum)
    call test_mpp_sum_5D(pe,npes,root,pelist,fullsum,pesum)
  if( pe.EQ.root ) print *, '------------------> Finished test_mpp_sum <------------------'

  call MPI_FINALIZE(ierr)

contains

  subroutine test_mpp_sum_scalar(pe,npes,root,pelist,fullsum,pesum)
    integer, intent(in) :: npes, pe, root, fullsum, pesum
    integer, intent(in), dimension(0:npes-2) :: pelist

    call test_mpp_sum_scalar_r4(pe,npes,root,pelist,fullsum,pesum)
    call test_mpp_sum_scalar_r8(pe,npes,root,pelist,fullsum,pesum)
    call test_mpp_sum_scalar_i4(pe,npes,root,pelist,fullsum,pesum)
    call test_mpp_sum_scalar_i8(pe,npes,root,pelist,fullsum,pesum)

  end subroutine test_mpp_sum_scalar

  !> Test the functionality of mpp_transmit for an r4_scalar.
  subroutine test_mpp_sum_scalar_r4(pe,npes,root,pelist,fullsum,pesum)
    integer, intent(in) :: npes, pe, root, fullsum, pesum
    integer, intent(in), dimension(0:npes-2) :: pelist
    real(kind=r4_kind) :: a4, suma4, b4, sumb4

    a4 = real(pe+1, kind=r4_kind)
    b4 = a4
    call mpp_sync()
    call mpp_sum(a4)
    if (pe .LE. npes-2) call mpp_sum(b4, pelist=pelist)
    if (pe .EQ. root) then
      suma4 = real(fullsum, kind=r4_kind)
      sumb4 = real(pesum, kind=r4_kind)
      if (a4 .ne. suma4) call mpp_error(FATAL, "Scalar_r4: mpp_sum differs from fortran intrinsic sum")
      if (b4 .ne. sumb4) call mpp_error(FATAL, "Scalar_r4 with pelist: mpp_sum differs from fortran intrinsic sum")
    endif

  end subroutine test_mpp_sum_scalar_r4

  !> Test the functionality of mpp_transmit for an r8_scalar.
  subroutine test_mpp_sum_scalar_r8(pe,npes,root,pelist,fullsum,pesum)
    integer, intent(in) :: npes, pe, root, fullsum, pesum
    integer, intent(in), dimension(0:npes-2) :: pelist
    real(kind=r8_kind) :: a8, suma8, b8, sumb8

    a8 = real(pe+1, kind=r8_kind)
    b8 = a8
    call mpp_sync()
    call mpp_sum(a8)
    if (pe .LE. npes-2) call mpp_sum(b8, pelist=pelist)
    if (pe .EQ. root) then
      suma8 = real(fullsum, kind=r8_kind)
      sumb8 = real(pesum, kind=r8_kind)
      if (a8 .ne. suma8) call mpp_error(FATAL, "Scalar_r8: mpp_sum differs from fortran intrinsic sum")
      if (b8 .ne. sumb8) call mpp_error(FATAL, "Scalar_r8 with pelist: mpp_sum differs from fortran intrinsic sum")
    endif

  end subroutine test_mpp_sum_scalar_r8

  !> Test the functionality of mpp_transmit for an i4_scalar.
  subroutine test_mpp_sum_scalar_i4(pe,npes,root,pelist,fullsum,pesum)
    integer, intent(in) :: npes, pe, root, fullsum, pesum
    integer, intent(in), dimension(0:npes-2) :: pelist
    integer(kind=i4_kind) :: a4, suma4, b4, sumb4

    a4 = int(pe+1, kind=i4_kind)
    b4 = a4
    call mpp_sync()
    call mpp_sum(a4)
    if (pe .LE. npes-2) call mpp_sum(b4, pelist=pelist)
    if (pe .EQ. root) then
      suma4 = int(fullsum, kind=i4_kind)
      sumb4 = int(pesum, kind=i4_kind)
      if (a4 .ne. suma4) call mpp_error(FATAL, "Scalar_i4: mpp_sum differs from fortran intrinsic sum")
      if (b4 .ne. sumb4) call mpp_error(FATAL, "Scalar_i4 with pelist: mpp_sum differs from fortran intrinsic sum")
    endif

  end subroutine test_mpp_sum_scalar_i4

  !> Test the functionality of mpp_transmit for an i8_scalar.
  subroutine test_mpp_sum_scalar_i8(pe,npes,root,pelist,fullsum,pesum)
    integer, intent(in) :: npes, pe, root, fullsum, pesum
    integer, intent(in), dimension(0:npes-2) :: pelist
    integer(kind=i8_kind) :: a8, suma8, b8, sumb8

    a8 = int(pe+1, kind=i8_kind)
    b8 = a8
    call mpp_sync()
    call mpp_sum(a8)
    if (pe .LE. npes-2) call mpp_sum(b8, pelist=pelist)
    if (pe .EQ. root) then
      suma8 = int(fullsum, kind=i8_kind)
      sumb8 = int(pesum, kind=i8_kind)
      if (a8 .ne. suma8) call mpp_error(FATAL, "Scalar_i8: mpp_sum differs from fortran intrinsic sum")
      if (b8 .ne. sumb8) call mpp_error(FATAL, "Scalar_i8 with pelist: mpp_sum differs from fortran intrinsic sum")
    endif

  end subroutine test_mpp_sum_scalar_i8

  subroutine test_mpp_sum_2D(pe,npes,root,pelist,fullsum,pesum)
    integer, intent(in) :: npes, pe, root, fullsum, pesum
    integer, intent(in), dimension(0:npes-2) :: pelist

    call test_mpp_sum_2D_r4(pe,npes,root,pelist,fullsum,pesum)
    call test_mpp_sum_2D_r8(pe,npes,root,pelist,fullsum,pesum)
    call test_mpp_sum_2D_i4(pe,npes,root,pelist,fullsum,pesum)
    call test_mpp_sum_2D_i8(pe,npes,root,pelist,fullsum,pesum)

  end subroutine test_mpp_sum_2D

  !> Test the functionality of mpp_transmit for an r4_2D.
  subroutine test_mpp_sum_2D_r4(pe,npes,root,pelist,fullsum,pesum)
    integer, intent(in) :: npes, pe, root, fullsum, pesum
    integer, intent(in), dimension(0:npes-2) :: pelist
    integer, parameter :: n=4
    real(kind=r4_kind), dimension(2,2) :: a4, suma4, b4, sumb4, c4, sumc4

    a4 = real(pe+1, kind=r4_kind)
    b4 = a4
    c4 = a4
    call mpp_sync()
    call mpp_sum(a4,n)
    if (pe .LE. npes-2) call mpp_sum(b4, n, pelist=pelist)
    call mpp_sum(c4,n-1)
    if (pe .EQ. root) then
      suma4 = real(fullsum, kind=r4_kind)
      sumb4 = real(pesum, kind=r4_kind)
      sumc4 = suma4
      sumc4(2,2) = 1
      if (all(a4 .ne. suma4)) call mpp_error(FATAL, "2D_r4: mpp_sum differs from fortran intrinsic sum")
      if (all(b4 .ne. sumb4)) call mpp_error(FATAL, "2D_r4 with pelist: mpp_sum differs from fortran intrinsic sum")
      if (all(c4 .ne. sumc4)) call mpp_error(FATAL, "2D_r4 (shorter length): mpp_sum differs from fortran intrinsic sum")
    endif

  end subroutine test_mpp_sum_2D_r4

  !> Test the functionality of mpp_transmit for an r8_2D.
  subroutine test_mpp_sum_2D_r8(pe,npes,root,pelist,fullsum,pesum)
    integer, intent(in) :: npes, pe, root, fullsum, pesum
    integer, intent(in), dimension(0:npes-2) :: pelist
    integer, parameter :: n=4
    real(kind=r8_kind), dimension(2,2) :: a8, suma8, b8, sumb8, c8, sumc8

    a8 = real(pe+1, kind=r8_kind)
    b8 = a8
    c8 = a8
    call mpp_sync()
    call mpp_sum(a8,n)
    if (pe .LE. npes-2) call mpp_sum(b8, n, pelist=pelist)
    call mpp_sum(c8, n-1)
    if (pe .EQ. root) then
      suma8 = real(fullsum, kind=r8_kind)
      sumb8 = real(pesum, kind=r8_kind)
      sumc8 = suma8
      sumc8(2,2) = 1
      if (all(a8 .ne. suma8)) call mpp_error(FATAL, "2D_r8: mpp_sum differs from fortran intrinsic sum")
      if (all(b8 .ne. sumb8)) call mpp_error(FATAL, "2D_r8 with pelist: mpp_sum differs from fortran intrinsic sum")
      if (all(c8 .ne. sumc8)) call mpp_error(FATAL, "2D_r8 (shorter length): mpp_sum differs from fortran intrinsic sum")
    endif

  end subroutine test_mpp_sum_2D_r8

  !> Test the functionality of mpp_transmit for an i4_2D.
  subroutine test_mpp_sum_2D_i4(pe,npes,root,pelist,fullsum,pesum)
    integer, intent(in) :: npes, pe, root, fullsum, pesum
    integer, intent(in), dimension(0:npes-2) :: pelist
    integer, parameter :: n=4
    integer(kind=i4_kind), dimension(2,2) :: a4, suma4, b4, sumb4, c4, sumc4

    a4 = int(pe+1, kind=i4_kind)
    b4 = a4
    c4 = a4
    call mpp_sync()
    call mpp_sum(a4,n)
    if (pe .LE. npes-2) call mpp_sum(b4, n, pelist=pelist)
    call mpp_sum(c4,n-1)
    if (pe .EQ. root) then
      suma4 = int(fullsum, kind=i4_kind)
      sumb4 = int(pesum, kind=i4_kind)
      sumc4 = suma4
      sumc4(2,2) = 1
      if (all(a4 .ne. suma4)) call mpp_error(FATAL, "2D_i4: mpp_sum differs from fortran intrinsic sum")
      if (all(b4 .ne. sumb4)) call mpp_error(FATAL, "2D_i4 with pelist: mpp_sum differs from fortran intrinsic sum")
      if (all(c4 .ne. sumc4)) call mpp_error(FATAL, "2D_i4 (shorter length): mpp_sum differs from fortran intrinsic sum")
    endif

  end subroutine test_mpp_sum_2D_i4

  !> Test the functionality of mpp_transmit for an i8_2D.
  subroutine test_mpp_sum_2D_i8(pe,npes,root,pelist,fullsum,pesum)
    integer, intent(in) :: npes, pe, root, fullsum, pesum
    integer, intent(in), dimension(0:npes-2) :: pelist
    integer, parameter :: n=4
    integer(kind=i8_kind), dimension(2,2) :: a8, suma8, b8, sumb8, c8, sumc8

    a8 = int(pe+1, kind=i8_kind)
    b8 = a8
    c8 = a8
    call mpp_sync()
    call mpp_sum(a8,n)
    if (pe .LE. npes-2) call mpp_sum(b8, n, pelist=pelist)
    call mpp_sum(c8,n-1)
    if (pe .EQ. root) then
      suma8 = int(fullsum, kind=i8_kind)
      sumb8 = int(pesum, kind=i8_kind)
      sumc8 = suma8
      sumc8(2,2) = 1
      if (all(a8 .ne. suma8)) call mpp_error(FATAL, "2D_i8: mpp_sum differs from fortran intrinsic sum")
      if (all(b8 .ne. sumb8)) call mpp_error(FATAL, "2D_i8 with pelist: mpp_sum differs from fortran intrinsic sum")
      if (all(c8 .ne. sumc8)) call mpp_error(FATAL, "2D_i8 (shorter length): mpp_sum differs from fortran intrinsic sum")
    endif

  end subroutine test_mpp_sum_2D_i8

  subroutine test_mpp_sum_3D(pe,npes,root,pelist,fullsum,pesum)
    integer, intent(in) :: npes, pe, root, fullsum, pesum
    integer, intent(in), dimension(0:npes-2) :: pelist

    call test_mpp_sum_3D_r4(pe,npes,root,pelist,fullsum,pesum)
    call test_mpp_sum_3D_r8(pe,npes,root,pelist,fullsum,pesum)
    call test_mpp_sum_3D_i4(pe,npes,root,pelist,fullsum,pesum)
    call test_mpp_sum_3D_i8(pe,npes,root,pelist,fullsum,pesum)

  end subroutine test_mpp_sum_3D

  !> Test the functionality of mpp_transmit for an r4_3D.
  subroutine test_mpp_sum_3D_r4(pe,npes,root,pelist,fullsum,pesum)
    integer, intent(in) :: npes, pe, root, fullsum, pesum
    integer, intent(in), dimension(0:npes-2) :: pelist
    integer, parameter :: n=8
    real(kind=r4_kind), dimension(2,2,2) :: a4, suma4, b4, sumb4, c4, sumc4

    a4 = real(pe+1, kind=r4_kind)
    b4 = a4
    c4 = a4
    call mpp_sync()
    call mpp_sum(a4,n)
    if (pe .LE. npes-2) call mpp_sum(b4, n, pelist=pelist)
    call mpp_sum(c4, n-1)
    if (pe .EQ. root) then
      suma4 = real(fullsum, kind=r4_kind)
      sumb4 = real(pesum, kind=r4_kind)
      sumc4 = suma4
      sumc4(2,2,2) = 1
      if (all(a4 .ne. suma4)) call mpp_error(FATAL, "3D_r4: mpp_sum differs from fortran intrinsic sum")
      if (all(b4 .ne. sumb4)) call mpp_error(FATAL, "3D_r4 with pelist: mpp_sum differs from fortran intrinsic sum")
      if (all(c4 .ne. sumc4)) call mpp_error(FATAL, "3D_r4 (shorter length): mpp_sum differs from fortran intrinsic sum")
    endif

  end subroutine test_mpp_sum_3D_r4

  !> Test the functionality of mpp_transmit for an r8_3D.
  subroutine test_mpp_sum_3D_r8(pe,npes,root,pelist,fullsum,pesum)
    integer, intent(in) :: npes, pe, root, fullsum, pesum
    integer, intent(in), dimension(0:npes-2) :: pelist
    integer, parameter :: n=8
    real(kind=r8_kind), dimension(2,2,2) :: a8, suma8, b8, sumb8, c8, sumc8

    a8 = real(pe+1, kind=r8_kind)
    b8 = a8
    c8 = a8
    call mpp_sync()
    call mpp_sum(a8,n)
    if (pe .LE. npes-2) call mpp_sum(b8, n, pelist=pelist)
    call mpp_sum(c8,n-1)
    if (pe .EQ. root) then
      suma8 = real(fullsum, kind=r8_kind)
      sumb8 = real(pesum, kind=r8_kind)
      sumc8 = suma8
      sumc8(2,2,2) = 1
      if (all(a8 .ne. suma8)) call mpp_error(FATAL, "3D_r8: mpp_sum differs from fortran intrinsic sum")
      if (all(b8 .ne. sumb8)) call mpp_error(FATAL, "3D_r8 with pelist: mpp_sum differs from fortran intrinsic sum")
      if (all(c8 .ne. sumc8)) call mpp_error(FATAL, "3D_r8 (shorter length): mpp_sum differs from fortran intrinsic sum")
    endif

  end subroutine test_mpp_sum_3D_r8

  !> Test the functionality of mpp_transmit for an i4_3D.
  subroutine test_mpp_sum_3D_i4(pe,npes,root,pelist,fullsum,pesum)
    integer, intent(in) :: npes, pe, root, fullsum, pesum
    integer, intent(in), dimension(0:npes-2) :: pelist
    integer, parameter :: n=8
    integer(kind=i4_kind), dimension(2,2,2) :: a4, suma4, b4, sumb4, c4, sumc4

    a4 = int(pe+1, kind=i4_kind)
    b4 = a4
    c4 = a4
    call mpp_sync()
    call mpp_sum(a4,n)
    if (pe .LE. npes-2) call mpp_sum(b4, n, pelist=pelist)
    call mpp_sum(c4,n-1)
    if (pe .EQ. root) then
      suma4 = int(fullsum, kind=i4_kind)
      sumb4 = int(pesum, kind=i4_kind)
      sumc4 = suma4
      sumc4(2,2,2) = 1
      if (all(a4 .ne. suma4)) call mpp_error(FATAL, "3D_i4: mpp_sum differs from fortran intrinsic sum")
      if (all(b4 .ne. sumb4)) call mpp_error(FATAL, "3D_i4 with pelist: mpp_sum differs from fortran intrinsic sum")
      if (all(c4 .ne. sumc4)) call mpp_error(FATAL, "3D_i4 (shorter length): mpp_sum differs from fortran intrinsic sum")
    endif

  end subroutine test_mpp_sum_3D_i4

  !> Test the functionality of mpp_transmit for an i8_3D.
  subroutine test_mpp_sum_3D_i8(pe,npes,root,pelist,fullsum,pesum)
    integer, intent(in) :: npes, pe, root, fullsum, pesum
    integer, intent(in), dimension(0:npes-2) :: pelist
    integer, parameter :: n=8
    integer(kind=i8_kind), dimension(2,2,2) :: a8, suma8, b8, sumb8, c8, sumc8

    a8 = int(pe+1, kind=i8_kind)
    b8 = a8
    c8 = a8
    call mpp_sync()
    call mpp_sum(a8,n)
    if (pe .LE. npes-2) call mpp_sum(b8, n, pelist=pelist)
    call mpp_sum(c8,n-1)
    if (pe .EQ. root) then
      suma8 = int(fullsum, kind=i8_kind)
      sumb8 = int(pesum, kind=i8_kind)
      sumc8 = suma8
      sumc8(2,2,2) = 1
      if (all(a8 .ne. suma8)) call mpp_error(FATAL, "3D_i8: mpp_sum differs from fortran intrinsic sum")
      if (all(b8 .ne. sumb8)) call mpp_error(FATAL, "3D_i8 with pelist: mpp_sum differs from fortran intrinsic sum")
      if (all(c8 .ne. sumc8)) call mpp_error(FATAL, "3D_i8 (shorter length): mpp_sum differs from fortran intrinsic sum")
    endif

  end subroutine test_mpp_sum_3D_i8

  subroutine test_mpp_sum_4D(pe,npes,root,pelist,fullsum,pesum)
    integer, intent(in) :: npes, pe, root, fullsum, pesum
    integer, intent(in), dimension(0:npes-2) :: pelist

    call test_mpp_sum_4D_r4(pe,npes,root,pelist,fullsum,pesum)
    call test_mpp_sum_4D_r8(pe,npes,root,pelist,fullsum,pesum)
    call test_mpp_sum_4D_i4(pe,npes,root,pelist,fullsum,pesum)
    call test_mpp_sum_4D_i8(pe,npes,root,pelist,fullsum,pesum)

  end subroutine test_mpp_sum_4D

  !> Test the functionality of mpp_transmit for an r4_4D.
  subroutine test_mpp_sum_4D_r4(pe,npes,root,pelist,fullsum,pesum)
    integer, intent(in) :: npes, pe, root, fullsum, pesum
    integer, intent(in), dimension(0:npes-2) :: pelist
    integer, parameter :: n=16
    real(kind=r4_kind), dimension(2,2,2,2) :: a4, suma4, b4, sumb4, c4, sumc4

    a4 = real(pe+1, kind=r4_kind)
    b4 = a4
    c4 = a4
    call mpp_sync()
    call mpp_sum(a4,n)
    if (pe .LE. npes-2) call mpp_sum(b4, n, pelist=pelist)
    call mpp_sum(c4,n-1)
    if (pe .EQ. root) then
      suma4 = real(fullsum, kind=r4_kind)
      sumb4 = real(pesum, kind=r4_kind)
      sumc4 = suma4
      sumc4(2,2,2,2) = 1
      if (all(a4 .ne. suma4)) call mpp_error(FATAL, "4D_r4: mpp_sum differs from fortran intrinsic sum")
      if (all(b4 .ne. sumb4)) call mpp_error(FATAL, "4D_r4 with pelist: mpp_sum differs from fortran intrinsic sum")
      if (all(c4 .ne. sumc4)) call mpp_error(FATAL, "4D_r4 (shorter length): mpp_sum differs from fortran intrinsic sum")
    endif

  end subroutine test_mpp_sum_4D_r4

  !> Test the functionality of mpp_transmit for an r8_4D.
  subroutine test_mpp_sum_4D_r8(pe,npes,root,pelist,fullsum,pesum)
    integer, intent(in) :: npes, pe, root, fullsum, pesum
    integer, intent(in), dimension(0:npes-2) :: pelist
    integer, parameter :: n=16
    real(kind=r8_kind), dimension(2,2,2,2) :: a8, suma8, b8, sumb8, c8, sumc8

    a8 = real(pe+1, kind=r8_kind)
    b8 = a8
    c8 = a8
    call mpp_sync()
    call mpp_sum(a8,n)
    if (pe .LE. npes-2) call mpp_sum(b8, n, pelist=pelist)
    call mpp_sum(c8,n-1)
    if (pe .EQ. root) then
      suma8 = real(fullsum, kind=r8_kind)
      sumb8 = real(pesum, kind=r8_kind)
      sumc8 = suma8
      sumc8(2,2,2,2) = 1
      if (all(a8 .ne. suma8)) call mpp_error(FATAL, "4D_r8: mpp_sum differs from fortran intrinsic sum")
      if (all(b8 .ne. sumb8)) call mpp_error(FATAL, "4D_r8 with pelist: mpp_sum differs from fortran intrinsic sum")
      if (all(c8 .ne. sumc8)) call mpp_error(FATAL, "4D_r8 (shorter length): mpp_sum differs from fortran intrinsic sum")
    endif

  end subroutine test_mpp_sum_4D_r8

  !> Test the functionality of mpp_transmit for an i4_4D.
  subroutine test_mpp_sum_4D_i4(pe,npes,root,pelist,fullsum,pesum)
    integer, intent(in) :: npes, pe, root, fullsum, pesum
    integer, intent(in), dimension(0:npes-2) :: pelist
    integer, parameter :: n=16
    integer(kind=i4_kind), dimension(2,2,2,2) :: a4, suma4, b4, sumb4, c4, sumc4

    a4 = int(pe+1, kind=i4_kind)
    b4 = a4
    c4 = a4
    call mpp_sync()
    call mpp_sum(a4,n)
    if (pe .LE. npes-2) call mpp_sum(b4, n, pelist=pelist)
    call mpp_sum(c4,n-1)
    if (pe .EQ. root) then
      suma4 = int(fullsum, kind=i4_kind)
      sumb4 = int(pesum, kind=i4_kind)
      sumc4 = suma4
      sumc4(2,2,2,2) = 1
      if (all(a4 .ne. suma4)) call mpp_error(FATAL, "4D_i4: mpp_sum differs from fortran intrinsic sum")
      if (all(b4 .ne. sumb4)) call mpp_error(FATAL, "4D_i4 with pelist: mpp_sum differs from fortran intrinsic sum")
      if (all(c4 .ne. sumc4)) call mpp_error(FATAL, "4D_i4 (shorter length): mpp_sum differs from fortran intrinsic sum")
    endif

  end subroutine test_mpp_sum_4D_i4

  !> Test the functionality of mpp_transmit for an i8_4D.
  subroutine test_mpp_sum_4D_i8(pe,npes,root,pelist,fullsum,pesum)
    integer, intent(in) :: npes, pe, root, fullsum, pesum
    integer, intent(in), dimension(0:npes-2) :: pelist
    integer, parameter :: n=16
    integer(kind=i8_kind), dimension(2,2,2,2) :: a8, suma8, b8, sumb8, c8, sumc8

    a8 = int(pe+1, kind=i8_kind)
    b8 = a8
    c8 = a8
    call mpp_sync()
    call mpp_sum(a8,n)
    if (pe .LE. npes-2) call mpp_sum(b8, n, pelist=pelist)
    call mpp_sum(c8,n-1)
    if (pe .EQ. root) then
      suma8 = int(fullsum, kind=i8_kind)
      sumb8 = int(pesum, kind=i8_kind)
      sumc8 = suma8
      sumc8(2,2,2,2) = 1
      if (all(a8 .ne. suma8)) call mpp_error(FATAL, "4D_i8: mpp_sum differs from fortran intrinsic sum")
      if (all(b8 .ne. sumb8)) call mpp_error(FATAL, "4D_i8 with pelist: mpp_sum differs from fortran intrinsic sum")
      if (all(c8 .ne. sumc8)) call mpp_error(FATAL, "4D_i8 (shorter length): mpp_sum differs from fortran intrinsic sum")
    endif

  end subroutine test_mpp_sum_4D_i8

  subroutine test_mpp_sum_5D(pe,npes,root,pelist,fullsum,pesum)
    integer, intent(in) :: npes, pe, root, fullsum, pesum
    integer, intent(in), dimension(0:npes-2) :: pelist

    call test_mpp_sum_5D_r4(pe,npes,root,pelist,fullsum,pesum)
    call test_mpp_sum_5D_r8(pe,npes,root,pelist,fullsum,pesum)
    call test_mpp_sum_5D_i4(pe,npes,root,pelist,fullsum,pesum)
    call test_mpp_sum_5D_i8(pe,npes,root,pelist,fullsum,pesum)

  end subroutine test_mpp_sum_5D

  !> Test the functionality of mpp_transmit for an r4_5D.
  subroutine test_mpp_sum_5D_r4(pe,npes,root,pelist,fullsum,pesum)
    integer, intent(in) :: npes, pe, root, fullsum, pesum
    integer, intent(in), dimension(0:npes-2) :: pelist
    integer, parameter :: n=32
    real(kind=r4_kind), dimension(2,2,2,2,2) :: a4, suma4, b4, sumb4, c4, sumc4

    a4 = real(pe+1, kind=r4_kind)
    b4 = a4
    c4 = a4
    call mpp_sync()
    call mpp_sum(a4,n)
    if (pe .LE. npes-2) call mpp_sum(b4, n, pelist=pelist)
    call mpp_sum(c4,n-1)
    if (pe .EQ. root) then
      suma4 = real(fullsum, kind=r4_kind)
      sumb4 = real(pesum, kind=r4_kind)
      sumc4 = suma4
      sumc4(2,2,2,2,2) = 1
      if (all(a4 .ne. suma4)) call mpp_error(FATAL, "5D_r4: mpp_sum differs from fortran intrinsic sum")
      if (all(b4 .ne. sumb4)) call mpp_error(FATAL, "5D_r4 with pelist: mpp_sum differs from fortran intrinsic sum")
      if (all(c4 .ne. sumc4)) call mpp_error(FATAL, "5D_r4 (shorter length): mpp_sum differs from fortran intrinsic sum")
    endif

  end subroutine test_mpp_sum_5D_r4

  !> Test the functionality of mpp_transmit for an r8_5D.
  subroutine test_mpp_sum_5D_r8(pe,npes,root,pelist,fullsum,pesum)
    integer, intent(in) :: npes, pe, root, fullsum, pesum
    integer, intent(in), dimension(0:npes-2) :: pelist
    integer, parameter :: n=32
    real(kind=r8_kind), dimension(2,2,2,2,2) :: a8, suma8, b8, sumb8, c8, sumc8

    a8 = real(pe+1, kind=r8_kind)
    b8 = a8
    c8 = a8
    call mpp_sync()
    call mpp_sum(a8,n)
    if (pe .LE. npes-2) call mpp_sum(b8, n, pelist=pelist)
    call mpp_sum(c8,n-1)
    if (pe .EQ. root) then
      suma8 = real(fullsum, kind=r8_kind)
      sumb8 = real(pesum, kind=r8_kind)
      sumc8 = suma8
      sumc8(2,2,2,2,2) = 1
      if (all(a8 .ne. suma8)) call mpp_error(FATAL, "5D_r8: mpp_sum differs from fortran intrinsic sum")
      if (all(b8 .ne. sumb8)) call mpp_error(FATAL, "5D_r8 with pelist: mpp_sum differs from fortran intrinsic sum")
      if (all(c8 .ne. sumc8)) call mpp_error(FATAL, "5D_r8 (shorter length): mpp_sum differs from fortran intrinsic sum")
    endif

  end subroutine test_mpp_sum_5D_r8

  !> Test the functionality of mpp_transmit for an i4_5D.
  subroutine test_mpp_sum_5D_i4(pe,npes,root,pelist,fullsum,pesum)
    integer, intent(in) :: npes, pe, root, fullsum, pesum
    integer, intent(in), dimension(0:npes-2) :: pelist
    integer, parameter :: n=32
    integer(kind=i4_kind), dimension(2,2,2,2,2) :: a4, suma4, b4, sumb4, c4, sumc4

    a4 = int(pe+1, kind=i4_kind)
    b4 = a4
    c4 = a4
    call mpp_sync()
    call mpp_sum(a4,n)
    if (pe .LE. npes-2) call mpp_sum(b4, n, pelist=pelist)
    call mpp_sum(c4,n-1)
    if (pe .EQ. root) then
      suma4 = int(fullsum, kind=i4_kind)
      sumb4 = int(pesum, kind=i4_kind)
      sumc4 = suma4
      sumc4(2,2,2,2,2) = 1
      if (all(a4 .ne. suma4)) call mpp_error(FATAL, "5D_i4: mpp_sum differs from fortran intrinsic sum")
      if (all(b4 .ne. sumb4)) call mpp_error(FATAL, "5D_i4 with pelist: mpp_sum differs from fortran intrinsic sum")
      if (all(c4 .ne. sumc4)) call mpp_error(FATAL, "5D_i4 (shorter length): mpp_sum differs from fortran intrinsic sum")
    endif

  end subroutine test_mpp_sum_5D_i4

  !> Test the functionality of mpp_transmit for an i8_5D.
  subroutine test_mpp_sum_5D_i8(pe,npes,root,pelist,fullsum,pesum)
    integer, intent(in) :: npes, pe, root, fullsum, pesum
    integer, intent(in), dimension(0:npes-2) :: pelist
    integer, parameter :: n=32
    integer(kind=i8_kind), dimension(2,2,2,2,2) :: a8, suma8, b8, sumb8, c8, sumc8

    a8 = int(pe+1, kind=i8_kind)
    b8 = a8
    c8 = a8
    call mpp_sync()
    call mpp_sum(a8,n)
    if (pe .LE. npes-2) call mpp_sum(b8, n, pelist=pelist)
    call mpp_sum(c8,n-1)
    if (pe .EQ. root) then
      suma8 = int(fullsum, kind=i8_kind)
      sumb8 = int(pesum, kind=i8_kind)
      sumc8 = suma8
      sumc8(2,2,2,2,2) = 1
      if (all(a8 .ne. suma8)) call mpp_error(FATAL, "5D_i8: mpp_sum differs from fortran intrinsic sum")
      if (all(b8 .ne. sumb8)) call mpp_error(FATAL, "5D_i8 with pelist: mpp_sum differs from fortran intrinsic sum")
      if (all(c8 .ne. sumc8)) call mpp_error(FATAL, "5D_i8 (shorter length): mpp_sum differs from fortran intrinsic sum")
    endif

  end subroutine test_mpp_sum_5D_i8
>>>>>>> b402a7097b2ec57cf3b0aafff80ccfad4773a20f

end program test_mpp_sum
