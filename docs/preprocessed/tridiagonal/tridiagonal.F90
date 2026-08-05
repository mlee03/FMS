# 1 "tridiagonal/tridiagonal.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "tridiagonal/tridiagonal.F90"
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
!> @defgroup tridiagonal_mod tridiagonal_mod
!! @ingroup tridiagonal
!! @{
!! @brief Solves a tridiagonal system of equations.
!!
!! The following schematic represents the system of equations solved,
!! where X is the solution.
!! <PRE>
!!     | B(1)  A(1)   0     0                .......            0    |  |X(1)|   |D(1)|
!!     | C(2)  B(2)  A(2)   0                .......            0    |  |X(2)|   |D(2)|
!!     |  0    C(3)  B(3)  A(3)  0           .......            0    |  | .. |   | .. |
!!     |  ..........................................                 |  | .. | = | .. |
!!     |  ..........................................                 |  | .. |   | .. |
!!     |                                  C(N-2) B(N-2) A(N-2)  0    |  | .. |   | .. |
!!     |                                    0    C(N-1) B(N-1) A(N-1)|  | .. |   | .. |
!!     |                                    0      0    C(N)   B(N)  |  |X(N)|   |D(N)|
!!
!! </PRE>
!!  To solve this system
!! <PRE>
!!   call tri_invert(X,D,A,B,C)
!!
!!       real, intent(out), dimension(:,:,:) :: X
!!       real, intent(in),  dimension(:,:,:) :: D
!!       real, optional,    dimension(:,:,:) :: A,B,C
!! </PRE>
!! For simplicity (?), A and C are assumed to be dimensioned the same size
!! as B, D, and X, although any input values for A(N) and C(1) are ignored.
!! (some checks are needed here)
!!
!! If A is not present, it is assumed that the matrix (A,B.C) has not been changed
!! since the last call to tri_invert.
!!
!! To release memory,
!! <PRE>
!!    call close_tridiagonal
!! </PRE>
!!
!!
!! Arguments A, B, and C are optional, and are saved as module variables
!! if one recalls tri_invert without changing (A,B,C)
!!
!! @note
!!     Optional arguments A,B,C have no intent declaration,
!!     so the default intent is inout. The value of A(N) is modified
!!     on output, and B and C are unchanged.
!!
!!  The following private allocatable arrays save the relevant information
!!  if one recalls tri_invert without changing (A,B,C):
!!  <PRE>
!!        allocate ( e  (size(x,1), size(x,2), size(x,3)) )
!!        allocate ( g  (size(x,1), size(x,2), size(x,3)) )
!!        allocate ( cc (size(x,1), size(x,2), size(x,3)) )
!!        allocate ( bb (size(x,1), size(x,2)) )
!! </PRE>
!!  This storage is deallocated when close_tridiagonal is called.

module tridiagonal_mod

    use platform_mod, only: r4_kind, r8_kind
    use mpp_mod,      only: mpp_error, FATAL
    implicit none

    type :: tridiag_reals_r4
        real(r4_kind), private, allocatable, dimension(:,:,:) :: e, g, cc
        real(r4_kind), private, allocatable, dimension(:,:)   :: bb
    end type

    type :: tridiag_reals_r8
        real(r8_kind), private, allocatable, dimension(:,:,:) :: e, g, cc
        real(r8_kind), private, allocatable, dimension(:,:)   :: bb
    end type

    type(tridiag_reals_r4) :: tridiag_r4 !< holds reals stored from r4_kind calls to tri_invert
    type(tridiag_reals_r8) :: tridiag_r8 !< holds reals stored from r8_kind calls to tri_invert

    logical, private :: init_tridiagonal_r4 = .false. !< true when fields in tridiag_r4 are allocated
    logical, private :: init_tridiagonal_r8 = .false. !< true when fields in tridiag_r8 are allocated

    !> Interface to solve tridiagonal systems of equations for either kind value.
    !! Module level variables will be deallocated and allocated for every
    !! Since this relies on the state of module variables (unless A,B,C are specified)
    !! the values stored are distinct for each kind call unless the added optional argument store_both_kinds
    !! is true.
    interface tri_invert
        module procedure tri_invert_r4
        module procedure tri_invert_r8
    end interface

    public :: tri_invert

    contains

    !> @brief Releases memory used by the solver
    subroutine close_tridiagonal
        if(.not. init_tridiagonal_r4 .and. .not. init_tridiagonal_r8) return
        !$OMP SINGLE
        if(allocated(tridiag_r4%e)) deallocate(tridiag_r4%e)
        if(allocated(tridiag_r4%g)) deallocate(tridiag_r4%g)
        if(allocated(tridiag_r4%cc)) deallocate(tridiag_r4%cc)
        if(allocated(tridiag_r4%bb)) deallocate(tridiag_r4%bb)
        if(allocated(tridiag_r8%e)) deallocate(tridiag_r8%e)
        if(allocated(tridiag_r8%g)) deallocate(tridiag_r8%g)
        if(allocated(tridiag_r8%cc)) deallocate(tridiag_r8%cc)
        if(allocated(tridiag_r8%bb)) deallocate(tridiag_r8%bb)
        init_tridiagonal_r4 = .false.; init_tridiagonal_r8 = .false.
        !$OMP END SINGLE
        return
    end subroutine close_tridiagonal


# 1 "tridiagonal/include/tridiagonal_r4.fh" 1
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************














# 1 "tridiagonal/include/tridiagonal.inc" 1
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************

!> @brief Sets up and solves the tridiagonal system of equations
!!
!! For simplicity, A and C are assumed to be dimensioned the same size
!! as B, D, and X, although any input values for A(N) and C(1) are ignored.
!! There are no checks to make sure the sizes agree.
!!
!! The value of A(N) is modified on output, and B and C are unchanged.
!!
!! For mixed precision, this routine uses the kind size macro(r4_kind) to determine
!! which module variables are used/stored. This means a,b, and c values will only be stored for calls
!! of the same real kind value unless store_both_kinds is present and .true..
subroutine tri_invert_r4(x,d,a,b,c, store_both_kinds)

    real(r4_kind), intent(out), dimension(:,:,:) :: x !< Solution to the tridiagonal system of equations
    real(r4_kind), intent(in),  dimension(:,:,:) :: d !< The right-hand side term, see the schematic above.
    real(r4_kind), optional,    dimension(:,:,:) :: a,b,c !< Left hand side terms(see schematic on module page).
                                                !! If not provided, values from last call are used
    logical, optional                                   :: store_both_kinds !< Will save module state
                                                         !! variables for both kind types in order to be used in
                                                         !! subsequent calls with either kind.

    real(r4_kind), dimension(size(x,1),size(x,2),size(x,3)) :: f
    integer, parameter :: kindl = r4_kind

    integer :: k

    if(present(a)) then
        !$OMP SINGLE
        init_tridiagonal_r4 = .true.
        if(allocated(tridiag_r4%e))     deallocate(tridiag_r4%e)
        if(allocated(tridiag_r4%g))     deallocate(tridiag_r4%g)
        if(allocated(tridiag_r4%bb))    deallocate(tridiag_r4%bb)
        if(allocated(tridiag_r4%cc))    deallocate(tridiag_r4%cc)
        allocate(tridiag_r4%e (size(x,1),size(x,2),size(x,3)))
        allocate(tridiag_r4%g (size(x,1),size(x,2),size(x,3)))
        allocate(tridiag_r4%bb(size(x,1),size(x,2)))
        allocate(tridiag_r4%cc(size(x,1),size(x,2),size(x,3)))
        !$OMP END SINGLE

        tridiag_r4%e(:,:,1) = - a(:,:,1) / b(:,:,1)
        a(:,:,size(x,3)) = 0.0_kindl

        do  k= 2,size(x,3)
            tridiag_r4%g(:,:,k) = 1.0_kindl/(b(:,:,k)+c(:,:,k)*tridiag_r4%e(:,:,k-1))
            tridiag_r4%e(:,:,k) = - a(:,:,k)* tridiag_r4%g(:,:,k)
        end do
        tridiag_r4%cc = c
        tridiag_r4%bb = 1.0_kindl/b(:,:,1)

    end if

    if(.not.init_tridiagonal_r4) call mpp_error(FATAL, 'tri_invert: a,b,and c args not provided or previously calculated.')

    f(:,:,1) =  d(:,:,1)*tridiag_r4%bb
    do k= 2, size(x,3)
        f(:,:,k) = (d(:,:,k) - tridiag_r4%cc(:,:,k)*f(:,:,k-1))*tridiag_r4%g(:,:,k)
    end do

    x(:,:,size(x,3)) = f(:,:,size(x,3))
    do k = size(x,3)-1,1,-1
        x(:,:,k) = tridiag_r4%e(:,:,k)*x(:,:,k+1)+f(:,:,k)
    end do

    ! stores both kind values for subsequent calculations if running with option
    if( present(store_both_kinds)) then
      if( store_both_kinds ) then
        if( r4_kind .eq. r8_kind) then
          tridiag_r4%e = real(tridiag_r4%e, r4_kind)
          tridiag_r4%g = real(tridiag_r4%g, r4_kind)
          tridiag_r4%cc = real(tridiag_r4%cc, r4_kind)
          tridiag_r4%bb = real(tridiag_r4%bb, r4_kind)
          init_tridiagonal_r4 = .true.
        else
          tridiag_r8%e = real(tridiag_r4%e, r8_kind)
          tridiag_r8%g = real(tridiag_r4%g, r8_kind)
          tridiag_r8%cc = real(tridiag_r4%cc, r8_kind)
          tridiag_r8%bb = real(tridiag_r4%bb, r8_kind)
          init_tridiagonal_r8 = .true.
        endif
      endif
    endif

    return
end subroutine tri_invert_r4
# 32 "tridiagonal/include/tridiagonal_r4.fh" 2
# 129 "tridiagonal/tridiagonal.F90" 2

# 1 "tridiagonal/include/tridiagonal_r8.fh" 1
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************














# 1 "tridiagonal/include/tridiagonal.inc" 1
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************

!> @brief Sets up and solves the tridiagonal system of equations
!!
!! For simplicity, A and C are assumed to be dimensioned the same size
!! as B, D, and X, although any input values for A(N) and C(1) are ignored.
!! There are no checks to make sure the sizes agree.
!!
!! The value of A(N) is modified on output, and B and C are unchanged.
!!
!! For mixed precision, this routine uses the kind size macro(r8_kind) to determine
!! which module variables are used/stored. This means a,b, and c values will only be stored for calls
!! of the same real kind value unless store_both_kinds is present and .true..
subroutine tri_invert_r8(x,d,a,b,c, store_both_kinds)

    real(r8_kind), intent(out), dimension(:,:,:) :: x !< Solution to the tridiagonal system of equations
    real(r8_kind), intent(in),  dimension(:,:,:) :: d !< The right-hand side term, see the schematic above.
    real(r8_kind), optional,    dimension(:,:,:) :: a,b,c !< Left hand side terms(see schematic on module page).
                                                !! If not provided, values from last call are used
    logical, optional                                   :: store_both_kinds !< Will save module state
                                                         !! variables for both kind types in order to be used in
                                                         !! subsequent calls with either kind.

    real(r8_kind), dimension(size(x,1),size(x,2),size(x,3)) :: f
    integer, parameter :: kindl = r8_kind

    integer :: k

    if(present(a)) then
        !$OMP SINGLE
        init_tridiagonal_r8 = .true.
        if(allocated(tridiag_r8%e))     deallocate(tridiag_r8%e)
        if(allocated(tridiag_r8%g))     deallocate(tridiag_r8%g)
        if(allocated(tridiag_r8%bb))    deallocate(tridiag_r8%bb)
        if(allocated(tridiag_r8%cc))    deallocate(tridiag_r8%cc)
        allocate(tridiag_r8%e (size(x,1),size(x,2),size(x,3)))
        allocate(tridiag_r8%g (size(x,1),size(x,2),size(x,3)))
        allocate(tridiag_r8%bb(size(x,1),size(x,2)))
        allocate(tridiag_r8%cc(size(x,1),size(x,2),size(x,3)))
        !$OMP END SINGLE

        tridiag_r8%e(:,:,1) = - a(:,:,1) / b(:,:,1)
        a(:,:,size(x,3)) = 0.0_kindl

        do  k= 2,size(x,3)
            tridiag_r8%g(:,:,k) = 1.0_kindl/(b(:,:,k)+c(:,:,k)*tridiag_r8%e(:,:,k-1))
            tridiag_r8%e(:,:,k) = - a(:,:,k)* tridiag_r8%g(:,:,k)
        end do
        tridiag_r8%cc = c
        tridiag_r8%bb = 1.0_kindl/b(:,:,1)

    end if

    if(.not.init_tridiagonal_r8) call mpp_error(FATAL, 'tri_invert: a,b,and c args not provided or previously calculated.')

    f(:,:,1) =  d(:,:,1)*tridiag_r8%bb
    do k= 2, size(x,3)
        f(:,:,k) = (d(:,:,k) - tridiag_r8%cc(:,:,k)*f(:,:,k-1))*tridiag_r8%g(:,:,k)
    end do

    x(:,:,size(x,3)) = f(:,:,size(x,3))
    do k = size(x,3)-1,1,-1
        x(:,:,k) = tridiag_r8%e(:,:,k)*x(:,:,k+1)+f(:,:,k)
    end do

    ! stores both kind values for subsequent calculations if running with option
    if( present(store_both_kinds)) then
      if( store_both_kinds ) then
        if( r8_kind .eq. r8_kind) then
          tridiag_r4%e = real(tridiag_r8%e, r4_kind)
          tridiag_r4%g = real(tridiag_r8%g, r4_kind)
          tridiag_r4%cc = real(tridiag_r8%cc, r4_kind)
          tridiag_r4%bb = real(tridiag_r8%bb, r4_kind)
          init_tridiagonal_r4 = .true.
        else
          tridiag_r8%e = real(tridiag_r8%e, r8_kind)
          tridiag_r8%g = real(tridiag_r8%g, r8_kind)
          tridiag_r8%cc = real(tridiag_r8%cc, r8_kind)
          tridiag_r8%bb = real(tridiag_r8%bb, r8_kind)
          init_tridiagonal_r8 = .true.
        endif
      endif
    endif

    return
end subroutine tri_invert_r8
# 32 "tridiagonal/include/tridiagonal_r8.fh" 2
# 130 "tridiagonal/tridiagonal.F90" 2

end module tridiagonal_mod

!> @}
! close documentation grouping
