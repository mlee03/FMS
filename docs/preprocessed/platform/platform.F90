# 1 "./platform/platform.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "./platform/platform.F90"
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
!> @defgroup platform_mod platform_mod
!! @ingroup platform
!! @{
!! @brief Uses @ref fms_platform.h to define byte sizes for variable kinds
!! to be used in fms.

module platform_mod
!platform-dependent settings

# 1 "./include/fms_platform.h" 1
! -*-f90-*-*
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





!Set type kinds.
# 37 "./include/fms_platform.h"
!These values are not necessarily portable.







!DEC$ MESSAGE:'Using 8-byte addressing'


!Control "pure" functions.





!DEC$ MESSAGE:'Using pure routines.'


!Control array members of derived types.
# 67 "./include/fms_platform.h"
!DEC$ MESSAGE:'Using allocatable derived type array members.'


!Control use of cray pointers within mpp_peset
!Other cray pointer usage in mpp routines is compiled regardless





!DEC$ MESSAGE:'Using cray pointers.'


 !If you do not want to use 64-bit integers.





!If you want to use quad-precision.





!Max string sizes for paths and files







# 27 "./platform/platform.F90" 2
  public
  integer, parameter :: r16_kind=8, r8_kind=8, r4_kind=4, &
                        c8_kind=8, c4_kind=4, &
                        l8_kind=8, l4_kind=4, &
                        i8_kind=8, i4_kind=4, i2_kind=2, &
                        ptr_kind=8
  integer, parameter :: FMS_PATH_LEN = 1024
  integer, parameter :: FMS_FILE_LEN = 255
!could additionally define things like OS, compiler...: useful?
end module platform_mod
!> @}
! close documentation grouping
