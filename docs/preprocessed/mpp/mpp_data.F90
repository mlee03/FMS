# 1 "mpp/mpp_data.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "mpp/mpp_data.F90"
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
!> @defgroup mpp_data_mod mpp_data_mod
!! @ingroup mpp
!! @{
!! @brief Module to hold pointer and stack data for use in @ref mpp modules.
!!
!! Makes stack and pointer data publicly available from @ref mpp_data_mpi.inc or @ref
!! mpp_data_nocomm.inc for use in @ref mpp modules. This module is mainly
!! for internal use within @ref mpp_mod and @ref mpp_domains_mod .

module mpp_data_mod





  use mpp_parameter_mod, only : MAXPES
  use platform_mod

  implicit none
  private

! Include variable "version" to be written to log file.

# 1 "./include/file_version.h" 1
! -*-f90-*-
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




  character(len=*), parameter :: version = 'unknown'
# 41 "mpp/mpp_data.F90" 2
  public version

  !> public data used by mpp_mod
  public :: stat, mpp_stack, ptr_stack, status, ptr_status, sync, ptr_sync
  public :: mpp_from_pe, ptr_from, remote_data_loc, ptr_remote

  !--- All othere modules should import these parameters from mpp_domains_mod.
  !> public data which is used by mpp_domains_mod.
  public :: mpp_domains_stack, ptr_domains_stack
  public :: mpp_domains_stack_nonblock, ptr_domains_stack_nonblock

  !-------------------------------------------------------------------------------!
  ! The following data included in the .inc file are diffrent for sma or mpi case !
  !-------------------------------------------------------------------------------!





# 1 "mpp/include/mpp_data_nocomm.inc" 1
! -*-f90-*-

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
!> @brief Holds dummy constants and stack data for @ref mpp_mod and @ref mpp_domains_mod.
!! Accessible through @ref mpp_data_mod

!----------------------------------------------------------------!
! The following data is used in mpp_mod and its components       !
!----------------------------------------------------------------!
real(r8_kind), allocatable :: mpp_stack(:)

!--- some dummy variables with dummy values that will never be used
integer, parameter :: stat=-999
integer, parameter :: ptr_stack = -999
integer, parameter :: status=-999, ptr_status = -999
integer, parameter :: remote_data_loc=-999, ptr_remote = -999
integer, parameter :: sync=-999, ptr_sync = -999
integer, parameter :: mpp_from_pe = -999, ptr_from = -999

!-------------------------------------------------------------------!
! The following data is used in mpp_domains_mod and its components  !
!-------------------------------------------------------------------!
real(r8_kind), allocatable, target :: mpp_domains_stack(:)
real(r8_kind), allocatable, target :: mpp_domains_stack_nonblock(:)
!--- some dummy variables with dummy values that will never be used
integer, parameter :: ptr_domains_stack = -999
integer, parameter :: ptr_domains_stack_nonblock = -999
# 60 "mpp/mpp_data.F90" 2


end module mpp_data_mod
!> @}
! close documentation grouping
