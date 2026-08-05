# 1 "horiz_interp/horiz_interp_type.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "horiz_interp/horiz_interp_type.F90"
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
!> @defgroup horiz_interp_type_mod horiz_interp_type_mod
!! @ingroup horiz_interp
!! @{
!! @brief define derived data type that contains indices and weights used for subsequent
!! interpolations.
!! @author Zhi Liang

module horiz_interp_type_mod

use mpp_mod, only : mpp_send, mpp_recv, mpp_sync_self, mpp_error, FATAL
use mpp_mod, only : mpp_pe, mpp_root_pe, mpp_npes
use mpp_mod, only : COMM_TAG_1, COMM_TAG_2
use platform_mod, only: r4_kind, r8_kind

implicit none
private


! parameter to determine interpolation method
 integer, parameter :: CONSERVE = 1
 integer, parameter :: BILINEAR = 2
 integer, parameter :: SPHERICAL = 3
 integer, parameter :: BICUBIC  = 4

public :: CONSERVE, BILINEAR, SPHERICAL, BICUBIC
public :: horiz_interp_type, stats, assignment(=)

interface assignment(=)
  module procedure horiz_interp_type_eq
end interface

interface stats
  module procedure stats_r4
  module procedure stats_r8
end interface


!> real(8) pointers for use in horiz_interp_type
type horizInterpReals8_type
   real(kind=r8_kind),    dimension(:,:), allocatable   :: faci     !< weights for conservative scheme
   real(kind=r8_kind),    dimension(:,:), allocatable   :: facj     !< weights for conservative scheme
   real(kind=r8_kind),    dimension(:,:), allocatable   :: area_src !< area of the source grid
   real(kind=r8_kind),    dimension(:,:), allocatable   :: area_dst !< area of the destination grid
   real(kind=r8_kind),    dimension(:,:,:), allocatable :: wti      !< weights for bilinear interpolation
                                                                    !! wti ist used for derivative "weights" in bicubic
   real(kind=r8_kind),    dimension(:,:,:), allocatable :: wtj      !< weights for bilinear interpolation
                                                                    !! wti ist used for derivative "weights" in bicubic
   real(kind=r8_kind),    dimension(:,:,:), allocatable :: src_dist !< distance between destination grid and
                                                                        !! neighbor source grid.
   real(kind=r8_kind),    dimension(:,:), allocatable   :: rat_x    !< the ratio of coordinates of the dest grid
                                                                    !! (x_dest -x_src_r)/(x_src_l -x_src_r)
                                                                    !! and (y_dest -y_src_r)/(y_src_l -y_src_r)
   real(kind=r8_kind),    dimension(:,:), allocatable   :: rat_y  !< the ratio of coordinates of the dest grid
                                                                  !! (x_dest -x_src_r)/(x_src_l -x_src_r)
                                                                  !! and (y_dest -y_src_r)/(y_src_l -y_src_r)
   real(kind=r8_kind),    dimension(:), allocatable     :: lon_in   !< the coordinates of the source grid
   real(kind=r8_kind),    dimension(:), allocatable     :: lat_in   !< the coordinates of the source grid
   real(kind=r8_kind),    dimension(:), allocatable     :: area_frac_dst !< area fraction in destination grid.
   real(kind=r8_kind),    dimension(:,:), allocatable   :: mask_in
   real(kind=r8_kind)                                   :: max_src_dist
   logical                                              :: is_allocated = .false. !< set to true upon field allocation

end type horizInterpReals8_type

!> holds real(4) pointers for use in horiz_interp_type
type horizInterpReals4_type
   real(kind=r4_kind),    dimension(:,:), allocatable   :: faci     !< weights for conservative scheme
   real(kind=r4_kind),    dimension(:,:), allocatable   :: facj     !< weights for conservative scheme
   real(kind=r4_kind),    dimension(:,:), allocatable   :: area_src !< area of the source grid
   real(kind=r4_kind),    dimension(:,:), allocatable   :: area_dst !< area of the destination grid
   real(kind=r4_kind),    dimension(:,:,:), allocatable :: wti      !< weights for bilinear interpolation
                                                                    !! wti ist used for derivative "weights" in bicubic
   real(kind=r4_kind),    dimension(:,:,:), allocatable :: wtj      !< weights for bilinear interpolation
                                                                    !! wti ist used for derivative "weights" in bicubic
   real(kind=r4_kind),    dimension(:,:,:), allocatable :: src_dist !< distance between destination grid and
                                                                        !! neighbor source grid.
   real(kind=r4_kind),    dimension(:,:), allocatable   :: rat_x    !< the ratio of coordinates of the dest grid
                                                                    !! (x_dest -x_src_r)/(x_src_l -x_src_r)
                                                                    !! and (y_dest -y_src_r)/(y_src_l -y_src_r)
   real(kind=r4_kind),    dimension(:,:), allocatable   :: rat_y  !< the ratio of coordinates of the dest grid
                                                                  !! (x_dest -x_src_r)/(x_src_l -x_src_r)
                                                                  !! and (y_dest -y_src_r)/(y_src_l -y_src_r)
   real(kind=r4_kind),    dimension(:), allocatable     :: lon_in   !< the coordinates of the source grid
   real(kind=r4_kind),    dimension(:), allocatable     :: lat_in   !< the coordinates of the source grid
   real(kind=r4_kind),    dimension(:), allocatable     :: area_frac_dst !< area fraction in destination grid.
   real(kind=r4_kind),    dimension(:,:), allocatable   :: mask_in
   real(kind=r4_kind)                                   :: max_src_dist
   logical                                              :: is_allocated = .false. !< set to true upon field allocation

end type horizInterpReals4_type

!> Holds data pointers and metadata for horizontal interpolations, passed between the horiz_interp modules
 type horiz_interp_type
   integer, dimension(:,:), allocatable   :: ilon    !< indices for conservative scheme
   integer, dimension(:,:), allocatable   :: jlat    !< indices for conservative scheme
                                                           !! wti ist used for derivative "weights" in bicubic
   integer, dimension(:,:,:), allocatable :: i_lon  !< indices for bilinear interpolation
                                                        !! and spherical regrid
   integer, dimension(:,:,:), allocatable :: j_lat  !< indices for bilinear interpolation
                                                        !! and spherical regrid
   logical, dimension(:,:), allocatable :: found_neighbors   !< indicate whether destination grid
                                                            !! has some source grid around it.
   integer, dimension(:,:), allocatable :: num_found
   integer                            :: nlon_src !< size of source grid
   integer                            :: nlat_src !< size of source grid
   integer                            :: nlon_dst !< size of destination grid
   integer                            :: nlat_dst !< size of destination grid
   integer                            :: interp_method      !< interpolation method.
                                                            !! =1, conservative scheme
                                                            !! =2, bilinear interpolation
                                                            !! =3, spherical regrid
                                                            !! =4, bicubic regrid
   logical                            :: I_am_initialized=.false.
   integer                            :: version                            !< indicate conservative
                                                                            !! interpolation version with value 1 or 2
   !--- The following are for conservative interpolation scheme version 2 ( through xgrid)
   integer                            :: nxgrid                             !< number of exchange grid
                                                                            !! between src and dst grid.
   integer, dimension(:), allocatable     :: i_src       !< indices in source grid.
   integer, dimension(:), allocatable     :: j_src       !< indices in source grid.
   integer, dimension(:), allocatable     :: i_dst       !< indices in destination grid.
   integer, dimension(:), allocatable     :: j_dst       !< indices in destination grid.
   type(horizInterpReals8_type) :: horizInterpReals8_type !< derived type holding kind 8 real data pointers
                                                                    !! if compiled with r8_kind
   type(horizInterpReals4_type) :: horizInterpReals4_type !< derived type holding kind 4 real data pointers
                                                                    !! if compiled with r8_kind
 end type

contains

!######################################################################################################################
!> @brief horiz_interp_type_eq creates a copy of the horiz_interp_type object
 subroutine horiz_interp_type_eq(horiz_interp_out, horiz_interp_in)
    type(horiz_interp_type), intent(inout) :: horiz_interp_out !< Output object being set
    type(horiz_interp_type), intent(in)    :: horiz_interp_in !< Input object being copied

    if(.not.horiz_interp_in%I_am_initialized) then
      call mpp_error(FATAL,'horiz_interp_type_eq: horiz_interp_type variable on right hand side is unassigned')
    endif

    if( allocated(horiz_interp_in%ilon )) &
      horiz_interp_out%ilon = horiz_interp_in%ilon

    if( allocated(horiz_interp_in%jlat )) &
      horiz_interp_out%jlat = horiz_interp_in%jlat

    if( allocated(horiz_interp_in%i_lon )) &
      horiz_interp_out%i_lon = horiz_interp_in%i_lon

    if( allocated(horiz_interp_in%j_lat )) &
      horiz_interp_out%j_lat = horiz_interp_in%j_lat

    if( allocated(horiz_interp_in%found_neighbors )) &
      horiz_interp_out%found_neighbors = horiz_interp_in%found_neighbors

    if( allocated(horiz_interp_in%num_found )) &
      horiz_interp_out%num_found = horiz_interp_in%num_found

    if( allocated(horiz_interp_in%i_src )) &
      horiz_interp_out%i_src = horiz_interp_in%i_src

    if( allocated(horiz_interp_in%j_src )) &
      horiz_interp_out%j_src = horiz_interp_in%j_src

    if( allocated(horiz_interp_in%i_dst )) &
      horiz_interp_out%i_dst = horiz_interp_in%i_dst

    if( allocated(horiz_interp_in%j_dst )) &
      horiz_interp_out%j_dst = horiz_interp_in%j_dst

    horiz_interp_out%nlon_src =  horiz_interp_in%nlon_src
    horiz_interp_out%nlat_src =  horiz_interp_in%nlat_src
    horiz_interp_out%nlon_dst =  horiz_interp_in%nlon_dst
    horiz_interp_out%nlat_dst =  horiz_interp_in%nlat_dst
    horiz_interp_out%interp_method   =  horiz_interp_in%interp_method
    horiz_interp_out%I_am_initialized = .true.

    if(horiz_interp_in%horizInterpReals8_type%is_allocated) then

      if( allocated(horiz_interp_in%horizInterpReals8_type%faci)) &
        horiz_interp_out%horizInterpReals8_type%faci = horiz_interp_in%horizInterpReals8_type%faci

      if( allocated(  horiz_interp_in%horizInterpReals8_type%facj)) &
        horiz_interp_out%horizInterpReals8_type%facj = horiz_interp_in%horizInterpReals8_type%facj

      if( allocated(  horiz_interp_in%horizInterpReals8_type%area_src)) &
        horiz_interp_out%horizInterpReals8_type%area_src = horiz_interp_in%horizInterpReals8_type%area_src

      if( allocated(  horiz_interp_in%horizInterpReals8_type%area_dst)) &
        horiz_interp_out%horizInterpReals8_type%area_dst = horiz_interp_in%horizInterpReals8_type%area_dst

      if( allocated(  horiz_interp_in%horizInterpReals8_type%wti)) &
        horiz_interp_out%horizInterpReals8_type%wti = horiz_interp_in%horizInterpReals8_type%wti

      if( allocated(  horiz_interp_in%horizInterpReals8_type%wtj)) &
        horiz_interp_out%horizInterpReals8_type%wtj = horiz_interp_in%horizInterpReals8_type%wtj

      if( allocated(  horiz_interp_in%horizInterpReals8_type%src_dist)) &
        horiz_interp_out%horizInterpReals8_type%src_dist = horiz_interp_in%horizInterpReals8_type%src_dist

      if( allocated(  horiz_interp_in%horizInterpReals8_type%rat_x)) &
        horiz_interp_out%horizInterpReals8_type%rat_x = horiz_interp_in%horizInterpReals8_type%rat_x

      if( allocated(  horiz_interp_in%horizInterpReals8_type%rat_y)) &
        horiz_interp_out%horizInterpReals8_type%rat_y = horiz_interp_in%horizInterpReals8_type%rat_y

      if( allocated(  horiz_interp_in%horizInterpReals8_type%lon_in)) &
        horiz_interp_out%horizInterpReals8_type%lon_in = horiz_interp_in%horizInterpReals8_type%lon_in

      if( allocated(  horiz_interp_in%horizInterpReals8_type%lat_in)) &
        horiz_interp_out%horizInterpReals8_type%lat_in = horiz_interp_in%horizInterpReals8_type%lat_in

      if( allocated(  horiz_interp_in%horizInterpReals8_type%area_frac_dst)) &
        horiz_interp_out%horizInterpReals8_type%area_frac_dst = horiz_interp_in%horizInterpReals8_type%area_frac_dst

      horiz_interp_out%horizInterpReals8_type%max_src_dist =  horiz_interp_in%horizInterpReals8_type%max_src_dist

      horiz_interp_out%horizInterpReals8_type%is_allocated = .true.
      ! this was left out previous to mixed mode
      if( allocated(horiz_interp_in%horizInterpReals8_type%mask_in)) &
        horiz_interp_out%horizInterpReals8_type%mask_in = horiz_interp_in%horizInterpReals8_type%mask_in

    else if (horiz_interp_in%horizInterpReals4_type%is_allocated) then
      if( allocated(horiz_interp_in%horizInterpReals4_type%faci)) &
        horiz_interp_out%horizInterpReals4_type%faci = horiz_interp_in%horizInterpReals4_type%faci

      if( allocated(  horiz_interp_in%horizInterpReals4_type%facj)) &
        horiz_interp_out%horizInterpReals4_type%facj = horiz_interp_in%horizInterpReals4_type%facj

      if( allocated(  horiz_interp_in%horizInterpReals4_type%area_src)) &
        horiz_interp_out%horizInterpReals4_type%area_src = horiz_interp_in%horizInterpReals4_type%area_src

      if( allocated(  horiz_interp_in%horizInterpReals4_type%area_dst)) &
        horiz_interp_out%horizInterpReals4_type%area_dst = horiz_interp_in%horizInterpReals4_type%area_dst

      if( allocated(  horiz_interp_in%horizInterpReals4_type%wti)) &
        horiz_interp_out%horizInterpReals4_type%wti = horiz_interp_in%horizInterpReals4_type%wti

      if( allocated(  horiz_interp_in%horizInterpReals4_type%wtj)) &
        horiz_interp_out%horizInterpReals4_type%wtj = horiz_interp_in%horizInterpReals4_type%wtj

      if( allocated(  horiz_interp_in%horizInterpReals4_type%src_dist)) &
        horiz_interp_out%horizInterpReals4_type%src_dist = horiz_interp_in%horizInterpReals4_type%src_dist

      if( allocated(  horiz_interp_in%horizInterpReals4_type%rat_x)) &
        horiz_interp_out%horizInterpReals4_type%rat_x = horiz_interp_in%horizInterpReals4_type%rat_x

      if( allocated(  horiz_interp_in%horizInterpReals4_type%rat_y)) &
        horiz_interp_out%horizInterpReals4_type%rat_y = horiz_interp_in%horizInterpReals4_type%rat_y

      if( allocated(  horiz_interp_in%horizInterpReals4_type%lon_in)) &
        horiz_interp_out%horizInterpReals4_type%lon_in = horiz_interp_in%horizInterpReals4_type%lon_in

      if( allocated(  horiz_interp_in%horizInterpReals4_type%lat_in)) &
        horiz_interp_out%horizInterpReals4_type%lat_in = horiz_interp_in%horizInterpReals4_type%lat_in

      if( allocated(  horiz_interp_in%horizInterpReals4_type%area_frac_dst)) &
        horiz_interp_out%horizInterpReals4_type%area_frac_dst  = horiz_interp_in%horizInterpReals4_type%area_frac_dst

      horiz_interp_out%horizInterpReals4_type%max_src_dist =  horiz_interp_in%horizInterpReals4_type%max_src_dist

      horiz_interp_out%horizInterpReals4_type%is_allocated = .true.
      ! this was left out previous to mixed mode
      if( allocated(horiz_interp_in%horizInterpReals4_type%mask_in)) &
        horiz_interp_out%horizInterpReals4_type%mask_in = horiz_interp_in%horizInterpReals4_type%mask_in

    else
        call mpp_error(FATAL, "horiz_interp_type_eq: cannot assign unallocated real values from horiz_interp_in")
    endif

    if(horiz_interp_in%interp_method == CONSERVE) then
        horiz_interp_out%version =  horiz_interp_in%version
        if(horiz_interp_in%version==2) horiz_interp_out%nxgrid = horiz_interp_in%nxgrid
    end if

 end subroutine horiz_interp_type_eq
!######################################################################################################################


# 1 "horiz_interp/include/horiz_interp_type_r4.fh" 1
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







# 1 "horiz_interp/include/horiz_interp_type.inc" 1
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
!> @brief This statistics is for bilinear interpolation and spherical regrid.
 subroutine stats_r4 ( dat, low, high, avg, miss, missing_value, mask )
 real(r4_kind),    intent(in)  :: dat(:,:)
 real(r4_kind),    intent(out) :: low, high, avg
 integer, intent(out) :: miss
 real(r4_kind), intent(in), optional :: missing_value
 real(r4_kind),    intent(in), optional :: mask(:,:)

 real(r4_kind) :: dsum, buffer_real(3)
 integer :: pe, root_pe, npes, p, buffer_int(2), npts
 integer, parameter :: kindl = r4_kind !< compiled kind size

   pe = mpp_pe()
   root_pe = mpp_root_pe()
   npes = mpp_npes()

   dsum = 0.0_kindl
   miss = 0

   if (present(missing_value)) then
      miss = count(dat(:,:) == missing_value)
      low  = minval(dat(:,:), dat(:,:) /= missing_value)
      high = maxval(dat(:,:), dat(:,:) /= missing_value)
      dsum = sum(dat(:,:), dat(:,:) /= missing_value)
   else if(present(mask)) then
      miss = count(mask(:,:) <= 0.5_kindl )
      low  = minval(dat(:,:),mask=mask(:,:) > 0.5_kindl)
      high = maxval(dat(:,:),mask=mask(:,:) > 0.5_kindl)
      dsum = sum(dat(:,:), mask=mask(:,:) > 0.5_kindl)
   else
      miss = 0
      low  = minval(dat(:,:))
      high = maxval(dat(:,:))
      dsum = sum(dat(:,:))
   endif
   avg = 0.0_kindl

   npts = size(dat(:,:)) - miss
   if(pe == root_pe) then
      do p = 1, npes - 1  ! root_pe receive data from other pe
      ! Force use of "scalar", integer pointer mpp interface
         call mpp_recv(buffer_real(1),glen=3, from_pe=p+root_pe, tag=COMM_TAG_1)
         dsum = dsum + buffer_real(1)
         low  = min(low, buffer_real(2))
         high = max(high, buffer_real(3))
         call mpp_recv(buffer_int(1), glen=2, from_pe=p+root_pe, tag=COMM_TAG_2)
         miss = miss + buffer_int(1)
         npts = npts + buffer_int(2)
      enddo
      if(npts == 0) then
         print*, 'Warning: no points is valid'
      else
         avg = dsum/real(npts, r4_kind)
      endif
    else   ! other pe send data to the root_pe.
      buffer_real(1) = dsum
      buffer_real(2) = low
      buffer_real(3) = high
      ! Force use of "scalar", integer pointer mpp interface
      call mpp_send(buffer_real(1),plen=3,to_pe=root_pe, tag=COMM_TAG_1)
      buffer_int(1) = miss
      buffer_int(2) = npts
      call mpp_send(buffer_int(1), plen=2, to_pe=root_pe, tag=COMM_TAG_2)
    endif

    call mpp_sync_self()

    return

 end subroutine stats_r4
# 25 "horiz_interp/include/horiz_interp_type_r4.fh" 2
# 297 "horiz_interp/horiz_interp_type.F90" 2

# 1 "horiz_interp/include/horiz_interp_type_r8.fh" 1
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







# 1 "horiz_interp/include/horiz_interp_type.inc" 1
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
!> @brief This statistics is for bilinear interpolation and spherical regrid.
 subroutine stats_r8 ( dat, low, high, avg, miss, missing_value, mask )
 real(r8_kind),    intent(in)  :: dat(:,:)
 real(r8_kind),    intent(out) :: low, high, avg
 integer, intent(out) :: miss
 real(r8_kind), intent(in), optional :: missing_value
 real(r8_kind),    intent(in), optional :: mask(:,:)

 real(r8_kind) :: dsum, buffer_real(3)
 integer :: pe, root_pe, npes, p, buffer_int(2), npts
 integer, parameter :: kindl = r8_kind !< compiled kind size

   pe = mpp_pe()
   root_pe = mpp_root_pe()
   npes = mpp_npes()

   dsum = 0.0_kindl
   miss = 0

   if (present(missing_value)) then
      miss = count(dat(:,:) == missing_value)
      low  = minval(dat(:,:), dat(:,:) /= missing_value)
      high = maxval(dat(:,:), dat(:,:) /= missing_value)
      dsum = sum(dat(:,:), dat(:,:) /= missing_value)
   else if(present(mask)) then
      miss = count(mask(:,:) <= 0.5_kindl )
      low  = minval(dat(:,:),mask=mask(:,:) > 0.5_kindl)
      high = maxval(dat(:,:),mask=mask(:,:) > 0.5_kindl)
      dsum = sum(dat(:,:), mask=mask(:,:) > 0.5_kindl)
   else
      miss = 0
      low  = minval(dat(:,:))
      high = maxval(dat(:,:))
      dsum = sum(dat(:,:))
   endif
   avg = 0.0_kindl

   npts = size(dat(:,:)) - miss
   if(pe == root_pe) then
      do p = 1, npes - 1  ! root_pe receive data from other pe
      ! Force use of "scalar", integer pointer mpp interface
         call mpp_recv(buffer_real(1),glen=3, from_pe=p+root_pe, tag=COMM_TAG_1)
         dsum = dsum + buffer_real(1)
         low  = min(low, buffer_real(2))
         high = max(high, buffer_real(3))
         call mpp_recv(buffer_int(1), glen=2, from_pe=p+root_pe, tag=COMM_TAG_2)
         miss = miss + buffer_int(1)
         npts = npts + buffer_int(2)
      enddo
      if(npts == 0) then
         print*, 'Warning: no points is valid'
      else
         avg = dsum/real(npts, r8_kind)
      endif
    else   ! other pe send data to the root_pe.
      buffer_real(1) = dsum
      buffer_real(2) = low
      buffer_real(3) = high
      ! Force use of "scalar", integer pointer mpp interface
      call mpp_send(buffer_real(1),plen=3,to_pe=root_pe, tag=COMM_TAG_1)
      buffer_int(1) = miss
      buffer_int(2) = npts
      call mpp_send(buffer_int(1), plen=2, to_pe=root_pe, tag=COMM_TAG_2)
    endif

    call mpp_sync_self()

    return

 end subroutine stats_r8
# 25 "horiz_interp/include/horiz_interp_type_r8.fh" 2
# 298 "horiz_interp/horiz_interp_type.F90" 2

end module horiz_interp_type_mod
!> @}
! close documentation grouping
