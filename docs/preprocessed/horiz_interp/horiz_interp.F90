# 1 "horiz_interp/horiz_interp.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "horiz_interp/horiz_interp.F90"
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
!> @defgroup horiz_interp_mod horiz_interp_mod
!! @ingroup horiz_interp
!! @{
!! @brief Performs spatial interpolation between grids.
!!
!! @author Zhi Liang, Bruce Wyman
!!
!! This module can interpolate data from any logically rectangular grid
!! to any logically rectangular grid. Four interpolation schems are used here:
!! conservative, bilinear, bicubic and inverse of square distance weighted.
!! The four interpolation schemes are implemented separately in
!! horiz_interp_conserver_mod, horiz_interp_blinear_mod, horiz_interp_bicubic_mod
!! and horiz_interp_spherical_mod. bicubic interpolation requires the source grid
!! is regular lon/lat grid. User can choose the interpolation method in the
!! public interface horiz_interp_new through optional argument interp_method,
!! with acceptable value "conservative", "bilinear", "bicubic" and "spherical".
!! The default value is "conservative". There is an optional mask field for
!! missing input data. An optional output mask field may be used in conjunction with
!! the input mask to show where output data exists.

module horiz_interp_mod

!-----------------------------------------------------------------------
!
!        Performs spatial interpolation between grids.
!
!-----------------------------------------------------------------------

use fms_mod,                    only: write_version_number, fms_error_handler
use fms_mod,                    only: check_nml_error
use mpp_mod,                    only: mpp_error, FATAL, stdout, stdlog, mpp_min
use mpp_mod,                    only: input_nml_file, WARNING, mpp_pe, mpp_root_pe
use constants_mod,              only: pi
use horiz_interp_type_mod,      only: horiz_interp_type, assignment(=)
use horiz_interp_type_mod,      only: CONSERVE, BILINEAR, SPHERICAL, BICUBIC
use horiz_interp_conserve_mod,  only: horiz_interp_conserve_init, horiz_interp_conserve
use horiz_interp_conserve_mod,  only: horiz_interp_conserve_new, horiz_interp_conserve_del
use horiz_interp_bilinear_mod,  only: horiz_interp_bilinear_init, horiz_interp_bilinear
use horiz_interp_bilinear_mod,  only: horiz_interp_bilinear_new, horiz_interp_bilinear_del
use horiz_interp_bilinear_mod,  only: horiz_interp_read_weights_bilinear
use horiz_interp_bicubic_mod,   only: horiz_interp_bicubic_init, horiz_interp_bicubic
use horiz_interp_bicubic_mod,   only: horiz_interp_bicubic_new, horiz_interp_bicubic_del
use horiz_interp_spherical_mod, only: horiz_interp_spherical_init, horiz_interp_spherical
use horiz_interp_spherical_mod, only: horiz_interp_spherical_new, horiz_interp_spherical_del
use platform_mod,               only: r4_kind, r8_kind

 implicit none
 private

!---- interfaces ----

 public   horiz_interp_type, horiz_interp, horiz_interp_new, horiz_interp_del, &
          horiz_interp_init, horiz_interp_end, assignment(=), horiz_interp_read_weights

!> Allocates space and initializes a derived-type variable
!! that contains pre-computed interpolation indices and weights.
!!
!> Allocates space and initializes a derived-type variable
!! that contains pre-computed interpolation indices and weights
!! for improved performance of multiple interpolations between
!! the same grids. This routine does not need to be called if you
!! are doing a single grid-to-grid interpolation.
!!
!! @param lon_in
!!      Longitude (in radians) for source data grid. You can pass 1-D lon_in to
!!      represent the geographical longitude of regular lon/lat grid, or just
!!      pass geographical longitude(lon_in is 2-D). The grid location may be
!!      located at grid cell edge or center, decided by optional argument "grid_at_center".
!!
!! @param lat_in
!!      Latitude (in radians) for source data grid. You can pass 1-D lat_in to
!!      represent the geographical latitude of regular lon/lat grid, or just
!!      pass geographical latitude(lat_in is 2-D). The grid location may be
!!      located at grid cell edge or center, decided by optional argument "grid_at_center".
!!
!! @param lon_out
!!      Longitude (in radians) for destination data grid. You can pass 1-D lon_out to
!!      represent the geographical longitude of regular lon/lat grid, or just
!!      pass geographical longitude(lon_out is 2-D). The grid location may be
!!      located at grid cell edge or center, decided by optional argument "grid_at_center".
!!
!! @param lat_out
!!      Latitude (in radians) for destination data grid. You can pass 1-D lat_out to
!!      represent the geographical latitude of regular lon/lat grid, or just
!!      pass geographical latitude(lat_out is 2-D). The grid location may be
!!      located at grid cell edge or center, decided by optional argument "grid_at_center".
!!
!! @param verbose
!!      Integer flag that controls the amount of printed output.
!!      verbose = 0, no output; = 1, min,max,means; = 2, still more
!!
!! @param interp_method
!!      interpolation method, = "conservative", using conservation scheme,
!!      = "bilinear", using bilinear interpolation, = "spherical",using spherical regrid.
!!      = "bicubic", using bicubic interpolation. The default value is "convervative".
!!
!! @param src_modulo
!!      Indicate the source data grid is cyclic or not.
!!
!! @param grid_at_center
!!      Indicate the data is on the center of grid box or the edge of grid box.
!!      When true, the data is on the center of grid box. default vaule is false.
!!      This option is only available when interp_method = "bilinear" or "bicubic".
!!
!! @param Interp
!!      A derived-type variable containing indices and weights used for subsequent
!!      interpolations. To reinitialize this variable for a different grid-to-grid
!!      interpolation you must first use the "horiz_interp_del" interface.
 interface horiz_interp_new
    ! Source grid is 1d, destination grid is 1d
    module procedure horiz_interp_new_1d_r4
    module procedure horiz_interp_new_1d_r8
    ! Source grid is 1d, destination grid is 2d
    module procedure horiz_interp_new_1d_src_r4
    module procedure horiz_interp_new_1d_src_r8
    ! Source grid is 2d, destination grid is 2d
    module procedure horiz_interp_new_2d_r4
    module procedure horiz_interp_new_2d_r8
    ! Source grid is 2d, destination grid is 1d
    module procedure horiz_interp_new_1d_dst_r4
    module procedure horiz_interp_new_1d_dst_r8
 end interface

 !> Subroutines for reading in weight files and using that to fill in the horiz_interp type instead
 !! calculating it
 interface horiz_interp_read_weights
   module procedure horiz_interp_read_weights_r4
   module procedure horiz_interp_read_weights_r8
 end interface horiz_interp_read_weights

!> Subroutine for performing the horizontal interpolation between two grids.
!!
!> Subroutine for performing the horizontal interpolation between
!! two grids. There are two forms of this interface.
!! Form A requires first calling horiz_interp_new, while Form B
!! requires no initialization.
!!
!! @param Interp
!!     Derived-type variable containing interpolation indices and weights.
!!     Returned by a previous call to horiz_interp_new.
!!
!! @param data_in
!!      Input data on source grid.
!!
!! @param verbose
!!      flag for the amount of print output.
!!               verbose = 0, no output; = 1, min,max,means; = 2, still more
!!
!! @param mask_in
!!      Input mask, must be the same size as the input data. The real value of
!!      mask_in must be in the range (0.,1.). Set mask_in=0.0 for data points
!!      that should not be used or have missing data. It is Not needed for
!!      spherical regrid.
!!
!! @param missing_value
!!      Use the missing_value to indicate missing data.
!!
!! @param missing_permit
!!      numbers of points allowed to miss for the bilinear interpolation. The value
!!      should be between 0 and 3.
!!
!! @param lon_in, lat_in
!!      longitude and latitude (in radians) of source grid. More explanation can
!!      be found in the documentation of horiz_interp_new.
!!
!! @param lon_out, lat_out
!!      longitude and latitude (in radians) of destination grid. More explanation can
!!      be found in the documentation of horiz_interp_new.
!!
!! @param data_out
!!      Output data on destination grid.
!!
!! @param mask_out
!!      Output mask that specifies whether data was computed.
!!
!!
!! @throws FATAL, size of input array incorrect
!!      The input data array does not match the size of the input grid edges
!!      specified. If you are using the initialization interface make sure you
!!      have the correct grid size.
!!
!! @throws FATAL, size of output array incorrect
!!      The output data array does not match the size of the input grid
!!      edges specified. If you are using the initialization interface make
!!      sure you have the correct grid size.
 interface horiz_interp
    module procedure horiz_interp_base_2d_r4
    module procedure horiz_interp_base_2d_r8
    module procedure horiz_interp_base_3d_r4
    module procedure horiz_interp_base_3d_r8
    module procedure horiz_interp_solo_1d_r4
    module procedure horiz_interp_solo_1d_r8
    module procedure horiz_interp_solo_1d_src_r4
    module procedure horiz_interp_solo_1d_src_r8
    module procedure horiz_interp_solo_2d_r4
    module procedure horiz_interp_solo_2d_r8
    module procedure horiz_interp_solo_1d_dst_r4
    module procedure horiz_interp_solo_1d_dst_r8
    module procedure horiz_interp_solo_old_r4
    module procedure horiz_interp_solo_old_r8
 end interface

!> Private helper routines
interface is_lat_lon
    module procedure is_lat_lon_r4
    module procedure is_lat_lon_r8
end interface

interface horiz_interp_solo_1d
  module procedure horiz_interp_solo_1d_r4
  module procedure horiz_interp_solo_1d_r8
end interface



 logical :: reproduce_siena = .false. !< Set reproduce_siena = .true. to reproduce siena results.
                 !! Set reproduce_siena = .false. to decrease truncation error
                 !! in routine poly_area in file mosaic_util.c. The truncation error of
                 !! second order conservative remapping might be big for high resolution
                 !! grid.

 namelist /horiz_interp_nml/ reproduce_siena

!-----------------------------------------------------------------------
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
# 244 "horiz_interp/horiz_interp.F90" 2
 logical            :: module_is_initialized = .FALSE.
!-----------------------------------------------------------------------

contains

!#######################################################################

  !> Initialize module and writes version number to logfile.out
  subroutine horiz_interp_init
  integer :: iunit, ierr, io

  if(module_is_initialized) return
  call write_version_number("HORIZ_INTERP_MOD", version)

  read (input_nml_file, horiz_interp_nml, iostat=io)
  ierr = check_nml_error(io,'horiz_interp_nml')
  if (mpp_pe() == mpp_root_pe() ) then
     iunit = stdlog()
     write (iunit, nml=horiz_interp_nml)
  endif

  if (reproduce_siena) then
    call mpp_error(FATAL, "horiz_interp_mod: You have overridden the default value of " // &
       "reproduce_siena and set it to .true. in horiz_interp_nml. This was a temporary workaround to " // &
       "allow for consistency in continuing experiments and is no longer supported. " // &
       "Please remove this namelist.")
  endif

  call horiz_interp_conserve_init
  call horiz_interp_bilinear_init
  call horiz_interp_bicubic_init
  call horiz_interp_spherical_init

  module_is_initialized = .true.

  end subroutine horiz_interp_init

!> Deallocates memory used by "horiz_interp_type" variables.
!! Must be called before reinitializing with horiz_interp_new.
 subroutine horiz_interp_del ( Interp )

   type (horiz_interp_type), intent(inout) :: Interp !< A derived-type variable returned by previous
                                           !! call to horiz_interp_new. The input variable must have
                                           !! allocated arrays. The returned variable will contain
                                           !! deallocated arrays

!-----------------------------------------------------------------------
!  releases space used by horiz_interp_type variables
!  must be called before re-initializing the same variable
!-----------------------------------------------------------------------
   select case(Interp % interp_method)
   case (CONSERVE)
      call horiz_interp_conserve_del(Interp )
   case (BILINEAR)
      call horiz_interp_bilinear_del(Interp )
   case (BICUBIC)
      call horiz_interp_bicubic_del(Interp )
   case (SPHERICAL)
      call horiz_interp_spherical_del(Interp )
   end select

   Interp%I_am_initialized = .false.
!-----------------------------------------------------------------------

 end subroutine horiz_interp_del

 !#####################################################################

 !> Dummy routine
 subroutine horiz_interp_end
 return
 end subroutine horiz_interp_end


# 1 "horiz_interp/include/horiz_interp_r4.fh" 1
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














































# 1 "horiz_interp/include/horiz_interp.inc" 1
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
  !> @brief Creates a 1D @ref horiz_interp_type with the given parameters
  subroutine horiz_interp_new_1d_r4 (Interp, lon_in, lat_in, lon_out, lat_out, verbose, &
                                  interp_method, num_nbrs, max_dist, src_modulo,     &
                                  grid_at_center, mask_in, mask_out)

    !-----------------------------------------------------------------------
    type(horiz_interp_type), intent(inout)        :: Interp
    real(r4_kind), intent(in),  dimension(:)               :: lon_in , lat_in
    real(r4_kind), intent(in),  dimension(:)               :: lon_out, lat_out
    integer, intent(in),                 optional :: verbose
    character(len=*), intent(in),        optional :: interp_method
    integer, intent(in),                 optional :: num_nbrs
    real(r4_kind),    intent(in),                 optional :: max_dist
    logical, intent(in),                 optional :: src_modulo
    logical, intent(in),                 optional :: grid_at_center
    real(r4_kind), intent(in), dimension(:,:),    optional :: mask_in  !< dummy variable
    real(r4_kind), intent(inout),dimension(:,:),  optional :: mask_out !< dummy variable
    !-----------------------------------------------------------------------
    real(r4_kind), dimension(:,:), allocatable :: lon_src, lat_src, lon_dst, lat_dst
    real(r4_kind), dimension(:),   allocatable :: lon_src_1d, lat_src_1d, lon_dst_1d, lat_dst_1d
    integer                           :: i, j, nlon_in, nlat_in, nlon_out, nlat_out
    logical                           :: center
    character(len=40)                 :: method
    integer, parameter                :: kindl = r4_kind !> real kind size currently compiling
    !-----------------------------------------------------------------------
    call horiz_interp_init

    method = 'conservative'
    if(present(interp_method)) method = interp_method

    select case (trim(method))
    case ("conservative")
       Interp%interp_method = CONSERVE
       call horiz_interp_conserve_new ( Interp, lon_in, lat_in, lon_out, lat_out, verbose)
    case ("bilinear")
       Interp%interp_method = BILINEAR
       center = .false.
       if(present(grid_at_center) ) center = grid_at_center
       if(center) then
          nlon_out = size(lon_out(:)); nlat_out = size(lat_out(:))
          allocate(lon_dst(nlon_out,nlat_out), lat_dst(nlon_out,nlat_out))
          do i = 1, nlon_out
             lon_dst(i,:) = lon_out(i)
          enddo
          do j = 1, nlat_out
             lat_dst(:,j) = lat_out(j)
          enddo

          call horiz_interp_bilinear_new ( Interp, lon_in, lat_in, lon_dst, lat_dst, &
               verbose, src_modulo)
          deallocate(lon_dst, lat_dst)
       else
          nlon_in  = size(lon_in(:))-1;  nlat_in  = size(lat_in(:))-1
          nlon_out = size(lon_out(:))-1; nlat_out = size(lat_out(:))-1
          allocate(lon_src_1d(nlon_in), lat_src_1d(nlat_in))
          allocate(lon_dst(nlon_out,nlat_out), lat_dst(nlon_out,nlat_out))
          do i = 1, nlon_in
             lon_src_1d(i) = (lon_in(i) + lon_in(i+1)) * 0.5_kindl
          enddo
          do j = 1, nlat_in
             lat_src_1d(j) = (lat_in(j) + lat_in(j+1)) * 0.5_kindl
          enddo
          do i = 1, nlon_out
             lon_dst(i,:) = (lon_out(i) + lon_out(i+1)) * 0.5_kindl
          enddo
          do j = 1, nlat_out
             lat_dst(:,j) = (lat_out(j) + lat_out(j+1)) * 0.5_kindl
          enddo
          call horiz_interp_bilinear_new ( Interp, lon_src_1d, lat_src_1d, lon_dst, lat_dst, &
               verbose, src_modulo)
          deallocate(lon_src_1d, lat_src_1d, lon_dst, lat_dst)
       endif
    case ("bicubic")
       Interp%interp_method = BICUBIC
       center = .false.
       if(present(grid_at_center) ) center = grid_at_center
       !No need to expand to 2d, horiz_interp_bicubic_new does 1d-1d
       if(center) then
          call horiz_interp_bicubic_new ( Interp, lon_in, lat_in, lon_out, lat_out, &
            verbose, src_modulo)
       else
          nlon_in  = size(lon_in(:))-1;  nlat_in  = size(lat_in(:))-1
          nlon_out = size(lon_out(:))-1; nlat_out = size(lat_out(:))-1
          allocate(lon_src_1d(nlon_in), lat_src_1d(nlat_in))
          allocate(lon_dst_1d(nlon_out), lat_dst_1d(nlat_out))
          do i = 1, nlon_in
             lon_src_1d(i) = (lon_in(i) + lon_in(i+1)) * 0.5_kindl
          enddo
          do j = 1, nlat_in
             lat_src_1d(j) = (lat_in(j) + lat_in(j+1)) * 0.5_kindl
          enddo
          do i = 1, nlon_out
             lon_dst_1d(i) = (lon_out(i) + lon_out(i+1)) * 0.5_kindl
          enddo
          do j = 1, nlat_out
             lat_dst_1d(j) = (lat_out(j) + lat_out(j+1)) * 0.5_kindl
          enddo
          call horiz_interp_bicubic_new ( Interp, lon_src_1d, lat_src_1d, lon_dst_1d, lat_dst_1d, &
               verbose, src_modulo)
          deallocate(lon_src_1d, lat_src_1d, lon_dst_1d, lat_dst_1d)
       endif
    case ("spherical")
       Interp%interp_method = SPHERICAL
       nlon_in  = size(lon_in(:));   nlat_in  = size(lat_in(:))
       nlon_out  = size(lon_out(:)); nlat_out = size(lat_out(:))
       allocate(lon_src(nlon_in,nlat_in), lat_src(nlon_in,nlat_in))
       allocate(lon_dst(nlon_out,nlat_out), lat_dst(nlon_out,nlat_out))
       do i = 1, nlon_in
          lon_src(i,:) = lon_in(i)
       enddo
       do j = 1, nlat_in
          lat_src(:,j) = lat_in(j)
       enddo
       do i = 1, nlon_out
          lon_dst(i,:) = lon_out(i)
       enddo
       do j = 1, nlat_out
          lat_dst(:,j) = lat_out(j)
       enddo
       call horiz_interp_spherical_new ( Interp, lon_src, lat_src, lon_dst, lat_dst, &
            num_nbrs, max_dist, src_modulo)
       deallocate(lon_src, lat_src, lon_dst, lat_dst)
    case default
       call mpp_error(FATAL,'horiz_interp_mod: interp_method should be conservative, bilinear, bicubic, spherical')
    end select

    !-----------------------------------------------------------------------
    Interp% horizInterpReals4_type % is_allocated = .true.
    Interp%I_am_initialized = .true.

  end subroutine horiz_interp_new_1d_r4

!#######################################################################

 subroutine horiz_interp_new_1d_src_r4 (Interp, lon_in, lat_in, lon_out, lat_out,   &
                                     verbose, interp_method, num_nbrs, max_dist, &
                                     src_modulo, grid_at_center, mask_in, mask_out, is_latlon_out )

   type(horiz_interp_type), intent(inout)        :: Interp
   real(r4_kind), intent(in),  dimension(:)               :: lon_in , lat_in
   real(r4_kind), intent(in),  dimension(:,:)             :: lon_out, lat_out
   integer, intent(in),                 optional :: verbose
   character(len=*), intent(in),        optional :: interp_method
   integer, intent(in),                 optional :: num_nbrs  !< minimum number of neighbors
   real(r4_kind),    intent(in),                 optional :: max_dist
   logical, intent(in),                 optional :: src_modulo
   logical, intent(in),                 optional :: grid_at_center
   real(r4_kind), intent(in), dimension(:,:),    optional :: mask_in
   real(r4_kind), intent(out),dimension(:,:),    optional :: mask_out
   logical, intent(in),                 optional :: is_latlon_out

   real(r4_kind), dimension(:,:), allocatable :: lon_src, lat_src
   real(r4_kind), dimension(:),   allocatable :: lon_src_1d, lat_src_1d
   integer                           :: i, j, nlon_in, nlat_in
   character(len=40)                 :: method
   logical                           :: center
   logical                           :: dst_is_latlon
   integer, parameter                :: kindl = r4_kind !< real kind size currently compiling
   !-----------------------------------------------------------------------
   call horiz_interp_init

   method = 'conservative'
   if(present(interp_method)) method = interp_method

   select case (trim(method))
   case ("conservative")
      Interp%interp_method = CONSERVE
      !--- check to see if the source grid is regular lat-lon grid or not.
      if(PRESENT(is_latlon_out)) then
         dst_is_latlon = is_latlon_out
      else
         dst_is_latlon = is_lat_lon(lon_out, lat_out)
      end if
      if(dst_is_latlon ) then
         if(present(mask_in)) then
            if ( ANY(mask_in < -.0001_kindl) .or. ANY(mask_in > 1.0001_kindl)) &
               call mpp_error(FATAL, &
                  'horiz_interp_conserve_new_1d_src(horiz_interp_conserve_mod): input mask not between 0,1')
            allocate(Interp%horizInterpReals4_type%mask_in(size(mask_in,1), size(mask_in,2)) )
            Interp%horizInterpReals4_type%mask_in = mask_in
         end if
         call horiz_interp_conserve_new ( Interp, lon_in, lat_in, lon_out(:,1), lat_out(1,:), &
              verbose=verbose )
      else
         call horiz_interp_conserve_new ( Interp, lon_in, lat_in, lon_out, lat_out, &
              verbose=verbose, mask_in=mask_in, mask_out=mask_out )
      end if
   case ("bilinear")
      Interp%interp_method = BILINEAR
      center = .false.
      if(present(grid_at_center) ) center = grid_at_center
      if(center) then
         call horiz_interp_bilinear_new ( Interp, lon_in, lat_in, lon_out, lat_out, &
              verbose, src_modulo )
      else
         nlon_in  = size(lon_in(:))-1;  nlat_in  = size(lat_in(:))-1
         allocate(lon_src_1d(nlon_in), lat_src_1d(nlat_in))
         do i = 1, nlon_in
            lon_src_1d(i) = (lon_in(i) + lon_in(i+1)) * 0.5_kindl
         enddo
         do j = 1, nlat_in
            lat_src_1d(j) = (lat_in(j) + lat_in(j+1)) * 0.5_kindl
         enddo
         call horiz_interp_bilinear_new ( Interp, lon_src_1d, lat_src_1d, lon_out, lat_out, &
              verbose, src_modulo )
         deallocate(lon_src_1d,lat_src_1d)
      endif
   case ("bicubic")
      Interp%interp_method = BICUBIC
      center = .false.
      if(present(grid_at_center) ) center = grid_at_center
      if(center) then
        call horiz_interp_bicubic_new ( Interp, lon_in, lat_in, lon_out, lat_out, &
              verbose, src_modulo )
      else
         nlon_in  = size(lon_in(:))-1;  nlat_in  = size(lat_in(:))-1
         allocate(lon_src_1d(nlon_in), lat_src_1d(nlat_in))
         do i = 1, nlon_in
            lon_src_1d(i) = (lon_in(i) + lon_in(i+1)) * 0.5_kindl
         enddo
         do j = 1, nlat_in
            lat_src_1d(j) = (lat_in(j) + lat_in(j+1)) * 0.5_kindl
         enddo
           call horiz_interp_bicubic_new ( Interp, lon_src_1d, lat_src_1d, lon_out, lat_out, &
              verbose, src_modulo )
         deallocate(lon_src_1d,lat_src_1d)
      endif
   case ("spherical")
      Interp%interp_method = SPHERICAL
      nlon_in  = size(lon_in(:));  nlat_in  = size(lat_in(:))
      allocate(lon_src(nlon_in,nlat_in), lat_src(nlon_in,nlat_in))
      do i = 1, nlon_in
         lon_src(i,:) = lon_in(i)
      enddo
      do j = 1, nlat_in
         lat_src(:,j) = lat_in(j)
      enddo
      call horiz_interp_spherical_new ( Interp, lon_src, lat_src, lon_out, lat_out, &
           num_nbrs, max_dist, src_modulo)
      deallocate(lon_src, lat_src)
   case default
      call mpp_error(FATAL,'interp_method should be conservative, bilinear, bicubic, spherical')
   end select

   !-----------------------------------------------------------------------
   Interp% horizInterpReals4_type % is_allocated = .true.
   Interp%I_am_initialized = .true.

 end subroutine horiz_interp_new_1d_src_r4

!#######################################################################

 subroutine horiz_interp_new_2d_r4 (Interp, lon_in, lat_in, lon_out, lat_out,   &
                                 verbose, interp_method, num_nbrs, max_dist, &
                                 src_modulo, mask_in, mask_out, is_latlon_in, is_latlon_out  )
 type(horiz_interp_type), intent(inout)     :: Interp
 real(r4_kind), intent(in),  dimension(:,:)          :: lon_in , lat_in
 real(r4_kind), intent(in),  dimension(:,:)          :: lon_out, lat_out
 integer, intent(in),              optional :: verbose
 character(len=*), intent(in),     optional :: interp_method
 integer, intent(in),              optional :: num_nbrs
 real(r4_kind),    intent(in),              optional :: max_dist
 logical, intent(in),              optional :: src_modulo
 real(r4_kind), intent(in), dimension(:,:), optional :: mask_in
 real(r4_kind), intent(out),dimension(:,:), optional :: mask_out
 logical, intent(in),              optional :: is_latlon_in, is_latlon_out
 logical           :: src_is_latlon, dst_is_latlon
 character(len=40) :: method
 integer, parameter                :: kindl = r4_kind !< real kind size currently compiling
!-----------------------------------------------------------------------
   call horiz_interp_init

   method = 'bilinear'
   if(present(interp_method)) method = interp_method

   select case (trim(method))
   case ("conservative")
      Interp%interp_method = CONSERVE
      if(PRESENT(is_latlon_in)) then
         src_is_latlon = is_latlon_in
      else
         src_is_latlon = is_lat_lon(lon_in, lat_in)
      end if
      if(PRESENT(is_latlon_out)) then
         dst_is_latlon = is_latlon_out
      else
         dst_is_latlon = is_lat_lon(lon_out, lat_out)
      end if
      if(src_is_latlon .AND. dst_is_latlon) then
         if(present(mask_in)) then
            if ( ANY(mask_in < -0.0001_kindl) .or. ANY(mask_in > 1.0001_kindl)) then
                call mpp_error(FATAL, 'horiz_interp_conserve_new_2d(horiz_interp_conserve_mod):' // &
                                       ' input mask not between 0,1')
            endif
            allocate(Interp%horizInterpReals4_type%mask_in(size(mask_in,1), size(mask_in,2)) )
            Interp%horizInterpReals4_type%mask_in = mask_in
         end if
         call horiz_interp_conserve_new ( Interp, lon_in(:,1), lat_in(1,:), lon_out(:,1), lat_out(1,:), &
              verbose=verbose )
      else if(src_is_latlon) then
         call horiz_interp_conserve_new ( Interp, lon_in(:,1), lat_in(1,:), lon_out, lat_out, &
              verbose=verbose, mask_in=mask_in, mask_out=mask_out )
      else if(dst_is_latlon) then
         call horiz_interp_conserve_new ( Interp, lon_in, lat_in, lon_out(:,1), lat_out(1,:), &
              verbose=verbose, mask_in=mask_in, mask_out=mask_out )
      else
         call horiz_interp_conserve_new ( Interp, lon_in, lat_in, lon_out, lat_out, &
              verbose=verbose, mask_in=mask_in, mask_out=mask_out )
      end if

   case ("spherical")
      Interp%interp_method = SPHERICAL
      call horiz_interp_spherical_new ( Interp, lon_in, lat_in, lon_out, lat_out, &
                                    num_nbrs, max_dist, src_modulo )
   case ("bilinear")
      Interp%interp_method = BILINEAR
      call horiz_interp_bilinear_new ( Interp, lon_in, lat_in, lon_out, lat_out, &
                                        verbose, src_modulo )
   case default
      call mpp_error(FATAL,'when source grid are 2d, interp_method should be spherical or bilinear')
   end select

   Interp% horizInterpReals4_type % is_allocated = .true.
   Interp%I_am_initialized = .true.

 end subroutine horiz_interp_new_2d_r4

!#######################################################################
 subroutine horiz_interp_new_1d_dst_r4 (Interp, lon_in, lat_in, lon_out, lat_out,   &
      verbose, interp_method, num_nbrs, max_dist, src_modulo, mask_in, mask_out, is_latlon_in )
   type(horiz_interp_type), intent(inout)     :: Interp
   real(r4_kind), intent(in),  dimension(:,:)          :: lon_in , lat_in
   real(r4_kind), intent(in),  dimension(:)            :: lon_out, lat_out
   integer, intent(in),              optional :: verbose
   character(len=*), intent(in),     optional :: interp_method
   integer, intent(in),              optional :: num_nbrs
   real(r4_kind),    intent(in),              optional :: max_dist
   logical, intent(in),              optional :: src_modulo
   real(r4_kind), intent(in), dimension(:,:), optional :: mask_in
   real(r4_kind), intent(out),dimension(:,:), optional :: mask_out
   logical, intent(in),              optional :: is_latlon_in

   character(len=40) :: method
   integer, parameter                :: kindl = r4_kind !< real kind size currently compiling
   !-------------some local variables-----------------------------------------------
   integer                           :: i, j, nlon_out, nlat_out
   real(r4_kind), dimension(:,:), allocatable :: lon_dst, lat_dst
   logical                           :: src_is_latlon
   !-----------------------------------------------------------------------
   call horiz_interp_init

   method = 'bilinear'
   if(present(interp_method)) method = interp_method

   nlon_out = size(lon_out(:)); nlat_out = size(lat_out(:))
   allocate(lon_dst(nlon_out,nlat_out), lat_dst(nlon_out,nlat_out))
   do i = 1, nlon_out
      lon_dst(i,:) = lon_out(i)
   enddo
   do j = 1, nlat_out
      lat_dst(:,j) = lat_out(j)
   enddo

   select case (trim(method))
   case ("conservative")
      Interp%interp_method = CONSERVE
      if(PRESENT(is_latlon_in)) then
         src_is_latlon = is_latlon_in
      else
         src_is_latlon = is_lat_lon(lon_in, lat_in)
      end if

      if(src_is_latlon) then
         if(present(mask_in)) then
            if ( ANY(mask_in < -0.0001_kindl) .or. ANY(mask_in > 1.0001_kindl)) &
                 call mpp_error(FATAL, &
                     'horiz_interp_conserve_new_1d_dst(horiz_interp_conserve_mod): input mask not between 0,1')
            allocate(Interp%horizInterpReals4_type%mask_in(size(mask_in,1), size(mask_in,2)) )
            Interp%horizInterpReals4_type%mask_in = mask_in
         end if
         call horiz_interp_conserve_new ( Interp, lon_in(:,1), lat_in(1,:), lon_out, lat_out, &
              verbose=verbose)
      else
         call horiz_interp_conserve_new ( Interp, lon_in, lat_in, lon_out, lat_out, &
              verbose=verbose, mask_in=mask_in, mask_out=mask_out )
      end if
   case ("bilinear")
      Interp%interp_method = BILINEAR
      call horiz_interp_bilinear_new ( Interp, lon_in, lat_in, lon_dst, lat_dst, &
           verbose, src_modulo )
   case ("spherical")
      Interp%interp_method = SPHERICAL
      call horiz_interp_spherical_new ( Interp, lon_in, lat_in, lon_dst, lat_dst, &
           num_nbrs, max_dist, src_modulo)
   case default
      call mpp_error(FATAL,'when source grid are 2d, interp_method should be spherical or bilinear')
   end select

   deallocate(lon_dst,lat_dst)

   Interp% horizInterpReals4_type % is_allocated = .true.
   Interp%I_am_initialized = .true.

 end subroutine horiz_interp_new_1d_dst_r4

!#######################################################################

 subroutine horiz_interp_base_2d_r4 ( Interp, data_in, data_out, verbose, &
                                   mask_in, mask_out, missing_value, missing_permit, &
                                   err_msg, new_missing_handle )
!-----------------------------------------------------------------------
   type (horiz_interp_type), intent(in) :: Interp
      real(r4_kind), intent(in),  dimension(:,:) :: data_in
      real(r4_kind), intent(out), dimension(:,:) :: data_out
   integer, intent(in),                   optional :: verbose
      real(r4_kind), intent(in),   dimension(:,:), optional :: mask_in
      real(r4_kind), intent(out),  dimension(:,:), optional :: mask_out
      real(r4_kind), intent(in),                   optional :: missing_value
      integer, intent(in),                optional :: missing_permit
   character(len=*), intent(out),         optional :: err_msg
      logical, intent(in),                optional :: new_missing_handle
!-----------------------------------------------------------------------
   if(present(err_msg)) err_msg = ''
   if(.not.Interp%I_am_initialized) then
     if(fms_error_handler('horiz_interp','The horiz_interp_type variable is not initialized',err_msg)) return
   endif

   select case(Interp%interp_method)
   case(CONSERVE)
      call horiz_interp_conserve(Interp,data_in, data_out, verbose, mask_in, mask_out)
   case(BILINEAR)
      call horiz_interp_bilinear(Interp,data_in, data_out, verbose, mask_in, mask_out, &
                             missing_value, missing_permit, new_missing_handle )
   case(BICUBIC)
      call horiz_interp_bicubic(Interp,data_in, data_out, verbose, mask_in, mask_out, &
                             missing_value, missing_permit )
   case(SPHERICAL)
      call horiz_interp_spherical(Interp,data_in, data_out, verbose, mask_in, mask_out, &
                             missing_value )
   case default
      call mpp_error(FATAL,'interp_method should be conservative, bilinear, bicubic, spherical')
   end select

   return

 end subroutine horiz_interp_base_2d_r4

!#######################################################################

 !> Overload of interface horiz_interp_base_2d_r4
 !! uses 3d arrays for data and mask
 !! this allows for multiple interpolations with one call
 subroutine horiz_interp_base_3d_r4 ( Interp, data_in, data_out, verbose, mask_in, mask_out, &
      missing_value, missing_permit, err_msg  )
   !-----------------------------------------------------------------------
   !   overload of interface horiz_interp_base_2d_r4
   !   uses 3d arrays for data and mask
   !   this allows for multiple interpolations with one call
   !-----------------------------------------------------------------------
   type (horiz_interp_type), intent(in)           :: Interp
   real(r4_kind), intent(in),  dimension(:,:,:)            :: data_in
   real(r4_kind), intent(out), dimension(:,:,:)            :: data_out
   integer, intent(in),                  optional :: verbose
   real(r4_kind), intent(in),   dimension(:,:,:), optional :: mask_in
   real(r4_kind), intent(out),  dimension(:,:,:), optional :: mask_out
   real(r4_kind), intent(in),                     optional :: missing_value
   integer, intent(in),                  optional :: missing_permit
   character(len=*), intent(out),        optional :: err_msg
   !-----------------------------------------------------------------------
   integer :: n

   if(present(err_msg)) err_msg = ''
   if(.not.Interp%I_am_initialized) then
     if(fms_error_handler('horiz_interp','The horiz_interp_type variable is not initialized',err_msg)) return
   endif

   do n = 1, size(data_in,3)
      if (present(mask_in))then
         if(present(mask_out)) then
            call horiz_interp( Interp, data_in(:,:,n), data_out(:,:,n), &
                 verbose, mask_in(:,:,n), mask_out(:,:,n), &
                 missing_value, missing_permit )
         else
            call horiz_interp( Interp, data_in(:,:,n), data_out(:,:,n), &
                 verbose, mask_in(:,:,n), missing_value = missing_value,  &
                 missing_permit = missing_permit )
         endif
      else
         if(present(mask_out)) then
            call horiz_interp( Interp, data_in(:,:,n), data_out(:,:,n), &
                 verbose, mask_out=mask_out(:,:,n), missing_value = missing_value,  &
                 missing_permit = missing_permit )
         else
            call horiz_interp( Interp, data_in(:,:,n), data_out(:,:,n), &
                 verbose, missing_value = missing_value,  &
                 missing_permit = missing_permit )
         endif
     endif
   enddo

   return
!-----------------------------------------------------------------------
 end subroutine horiz_interp_base_3d_r4

!#######################################################################

!> Interpolates from a rectangular grid to rectangular grid.
!! interp_method can be the value conservative, bilinear or spherical.
!! horiz_interp_new don't need to be called before calling this routine.
 subroutine horiz_interp_solo_1d_r4 ( data_in, lon_in, lat_in, lon_out, lat_out,    &
                                   data_out, verbose, mask_in, mask_out,         &
                                   interp_method, missing_value, missing_permit, &
                                   num_nbrs, max_dist,src_modulo, grid_at_center  )
!-----------------------------------------------------------------------
      real(r4_kind), intent(in),  dimension(:,:) :: data_in
      real(r4_kind), intent(in),  dimension(:)   :: lon_in , lat_in
      real(r4_kind), intent(in),  dimension(:)   :: lon_out, lat_out
      real(r4_kind), intent(out), dimension(:,:) :: data_out
   integer, intent(in),                   optional :: verbose
      real(r4_kind), intent(in),   dimension(:,:), optional :: mask_in
      real(r4_kind), intent(out),  dimension(:,:), optional :: mask_out
   character(len=*), intent(in),          optional :: interp_method
      real(r4_kind), intent(in),                   optional :: missing_value
   integer, intent(in),                   optional :: missing_permit
   integer, intent(in),                   optional :: num_nbrs
      real(r4_kind), intent(in),                   optional :: max_dist
   logical, intent(in),                   optional :: src_modulo
   logical, intent(in),                   optional :: grid_at_center
!-----------------------------------------------------------------------
    type (horiz_interp_type) :: Interp
!-----------------------------------------------------------------------
    call horiz_interp_init

    call horiz_interp_new ( Interp, lon_in, lat_in, lon_out, lat_out, verbose, &
                             interp_method, num_nbrs, max_dist, src_modulo, grid_at_center )

    call horiz_interp ( Interp, data_in, data_out, verbose,   &
                        mask_in, mask_out, missing_value, missing_permit )

    call horiz_interp_del ( Interp )
!-----------------------------------------------------------------------

 end subroutine horiz_interp_solo_1d_r4

!#######################################################################

!> Interpolates from a uniformly spaced grid to any output grid.
!! interp_method can be the value "onservative","bilinear" or "spherical".
!! horiz_interp_new don't need to be called before calling this routine.
 subroutine horiz_interp_solo_1d_src_r4 ( data_in, lon_in, lat_in, lon_out, lat_out,    &
                                       data_out, verbose, mask_in, mask_out,         &
                                       interp_method, missing_value, missing_permit, &
                                       num_nbrs, max_dist, src_modulo, grid_at_center )
!-----------------------------------------------------------------------
      real(r4_kind), intent(in),  dimension(:,:) :: data_in
      real(r4_kind), intent(in),  dimension(:)   :: lon_in , lat_in
      real(r4_kind), intent(in),  dimension(:,:) :: lon_out, lat_out
      real(r4_kind), intent(out), dimension(:,:) :: data_out
   integer, intent(in),                   optional :: verbose
      real(r4_kind), intent(in),   dimension(:,:), optional :: mask_in
      real(r4_kind), intent(out),  dimension(:,:), optional :: mask_out
   character(len=*), intent(in),          optional :: interp_method
      real(r4_kind), intent(in),                   optional :: missing_value
   integer, intent(in),                   optional :: missing_permit
   integer, intent(in),                   optional :: num_nbrs
      real(r4_kind), intent(in),                   optional :: max_dist
   logical, intent(in),                   optional :: src_modulo
   logical, intent(in),                   optional :: grid_at_center

!-----------------------------------------------------------------------
   type (horiz_interp_type) :: Interp
   logical                  :: dst_is_latlon
   character(len=128)       :: method
!-----------------------------------------------------------------------
    call horiz_interp_init
    method = 'conservative'
    if(present(interp_method)) method = interp_method
    dst_is_latlon = .true.
    if(trim(method) == 'conservative') dst_is_latlon = is_lat_lon(lon_out, lat_out)

    if(dst_is_latlon) then
       call horiz_interp_new ( Interp, lon_in, lat_in, lon_out, lat_out, verbose, &
                               interp_method, num_nbrs, max_dist, src_modulo,    &
                               grid_at_center, is_latlon_out = dst_is_latlon )
       call horiz_interp ( Interp, data_in, data_out, verbose,   &
                           mask_in, mask_out, missing_value, missing_permit )
    else
       call horiz_interp_new ( Interp, lon_in, lat_in, lon_out, lat_out, verbose, &
                               interp_method, num_nbrs, max_dist, src_modulo,    &
                               grid_at_center, mask_in, mask_out, is_latlon_out = dst_is_latlon)

       call horiz_interp ( Interp, data_in, data_out, verbose,   &
                           missing_value=missing_value, missing_permit=missing_permit )
    end if

    call horiz_interp_del ( Interp )

!-----------------------------------------------------------------------

 end subroutine horiz_interp_solo_1d_src_r4


!#######################################################################

!> Interpolates from any grid to any grid. interp_method should be "spherical"
!! horiz_interp_new don't need to be called before calling this routine.
 subroutine horiz_interp_solo_2d_r4 ( data_in, lon_in, lat_in, lon_out, lat_out, data_out, &
                                   verbose, mask_in, mask_out, interp_method, missing_value,&
                                   missing_permit, num_nbrs, max_dist, src_modulo  )
!-----------------------------------------------------------------------
      real(r4_kind), intent(in),  dimension(:,:) :: data_in
      real(r4_kind), intent(in),  dimension(:,:) :: lon_in , lat_in
      real(r4_kind), intent(in),  dimension(:,:) :: lon_out, lat_out
      real(r4_kind), intent(out), dimension(:,:) :: data_out
   integer, intent(in),                   optional :: verbose
      real(r4_kind), intent(in),   dimension(:,:), optional :: mask_in
      real(r4_kind), intent(out),  dimension(:,:), optional :: mask_out
   character(len=*), intent(in),          optional :: interp_method
      real(r4_kind), intent(in),                   optional :: missing_value
   integer, intent(in),                   optional :: missing_permit
   integer, intent(in),                   optional :: num_nbrs
      real(r4_kind), intent(in),                   optional :: max_dist
   logical, intent(in),                   optional :: src_modulo
!-----------------------------------------------------------------------
   type (horiz_interp_type) :: Interp
   logical                  :: dst_is_latlon, src_is_latlon
   character(len=128)       :: method
!-----------------------------------------------------------------------
    call horiz_interp_init

    method = 'conservative'
    if(present(interp_method)) method = interp_method
    dst_is_latlon = .true.
    src_is_latlon = .true.
    if(trim(method) == 'conservative') then
       dst_is_latlon = is_lat_lon(lon_out, lat_out)
       src_is_latlon = is_lat_lon(lon_in, lat_in)
    end if

    if(dst_is_latlon .and. src_is_latlon) then
       call horiz_interp_new ( Interp, lon_in, lat_in, lon_out, lat_out, verbose, &
                               interp_method, num_nbrs, max_dist, src_modulo,    &
                               is_latlon_in=dst_is_latlon, is_latlon_out = dst_is_latlon )
       call horiz_interp ( Interp, data_in, data_out, verbose,   &
                           mask_in, mask_out, missing_value, missing_permit )
    else
       call horiz_interp_new ( Interp, lon_in, lat_in, lon_out, lat_out, verbose, &
                               interp_method, num_nbrs, max_dist, src_modulo,    &
                               mask_in, mask_out, &
                               is_latlon_in=dst_is_latlon, is_latlon_out = dst_is_latlon)
       call horiz_interp ( Interp, data_in, data_out, verbose,   &
                           missing_value=missing_value, missing_permit=missing_permit )
    end if

    call horiz_interp_del ( Interp )

!-----------------------------------------------------------------------

 end subroutine horiz_interp_solo_2d_r4

!#######################################################################

!>   interpolates from any grid to rectangular longitude/latitude grid.
!!   interp_method should be "spherical".
!!   horiz_interp_new don't need to be called before calling this routine.
 subroutine horiz_interp_solo_1d_dst_r4 ( data_in, lon_in, lat_in, lon_out, lat_out, data_out,    &
                                       verbose, mask_in, mask_out,interp_method,missing_value, &
                                       missing_permit,  num_nbrs, max_dist, src_modulo)
!-----------------------------------------------------------------------
      real(r4_kind), intent(in),  dimension(:,:) :: data_in
      real(r4_kind), intent(in),  dimension(:,:) :: lon_in , lat_in
      real(r4_kind), intent(in),  dimension(:)   :: lon_out, lat_out
      real(r4_kind), intent(out), dimension(:,:) :: data_out
   integer, intent(in),                   optional :: verbose
      real(r4_kind), intent(in),   dimension(:,:), optional :: mask_in
      real(r4_kind), intent(out),  dimension(:,:), optional :: mask_out
   character(len=*), intent(in),          optional :: interp_method
      real(r4_kind), intent(in),                   optional :: missing_value
   integer, intent(in),                   optional :: missing_permit
   integer, intent(in),                   optional :: num_nbrs
      real(r4_kind), intent(in),                   optional :: max_dist
   logical, intent(in),                   optional :: src_modulo
!-----------------------------------------------------------------------
   type (horiz_interp_type) :: Interp
   logical                  :: src_is_latlon
   character(len=128)       :: method
!-----------------------------------------------------------------------
    call horiz_interp_init

    method = 'conservative'
    if(present(interp_method)) method = interp_method
    src_is_latlon = .true.
    if(trim(method) == 'conservative') src_is_latlon = is_lat_lon(lon_in, lat_in)

    if(src_is_latlon) then
       call horiz_interp_new ( Interp, lon_in, lat_in, lon_out, lat_out, verbose, &
                               interp_method, num_nbrs, max_dist, src_modulo,    &
                               is_latlon_in = src_is_latlon )
       call horiz_interp ( Interp, data_in, data_out, verbose,   &
                           mask_in, mask_out, missing_value, missing_permit )
    else
       call horiz_interp_new ( Interp, lon_in, lat_in, lon_out, lat_out, verbose, &
                               interp_method, num_nbrs, max_dist, src_modulo,    &
                               mask_in, mask_out, is_latlon_in = src_is_latlon)

       call horiz_interp ( Interp, data_in, data_out, verbose,   &
                           missing_value=missing_value, missing_permit=missing_permit )
    end if

    call horiz_interp_del ( Interp )

!-----------------------------------------------------------------------

 end subroutine horiz_interp_solo_1d_dst_r4

!#######################################################################

!> Overloaded version of interface horiz_interp_solo_2
 subroutine horiz_interp_solo_old_r4 (data_in, wb, sb, dx, dy,  &
                                   lon_out, lat_out, data_out,  &
                                   verbose, mask_in, mask_out)

!-----------------------------------------------------------------------
      real(r4_kind), intent(in),  dimension(:,:) :: data_in !< Global input data stored from west to east
                                        !! (1st dimension), south to north (2nd dimension)
      real(r4_kind), intent(in)                  :: wb !< Longitude (radians) that correspond to western-most
                                              !! boundary of grid box j=1 in array data_in
      real(r4_kind), intent(in)                  :: sb !< Latitude (radians) that correspond to western-most
                                              !! boundary of grid box j=1 in array data_in
      real(r4_kind), intent(in)                  :: dx !< Grid spacing (in radians) for the longitude axis
                                              !! (first dimension) for the input data
      real(r4_kind), intent(in)                  :: dy !< Grid spacing (in radians) for the latitude axis
                                              !! (first dimension) for the input data
      real(r4_kind), intent(in),  dimension(:)   :: lon_out !< The longitude edges (in radians) for output
                                        !! data grid boxes. The values are for adjacent grid boxes
                                        !! and must increase in value. If there are MLON grid boxes
                                        !! there must be MLON+1 edge values
      real(r4_kind), intent(in),  dimension(:)   :: lat_out !< The latitude edges (in radians) for output
                                        !! data grid boxes. The values are for adjacent grid boxes
                                        !! and may increase or decrease in value. If there are NLAT
                                        !! grid boxes there must be NLAT+1 edge values
      real(r4_kind), intent(out), dimension(:,:) :: data_out !< Output data on the output grid defined by grid box
   integer, intent(in),                   optional :: verbose
      real(r4_kind), intent(in),   dimension(:,:), optional :: mask_in
      real(r4_kind), intent(out),  dimension(:,:), optional :: mask_out
!-----------------------------------------------------------------------
     real(r4_kind), dimension(size(data_in,1)+1)  :: blon_in
     real(r4_kind), dimension(size(data_in,2)+1)  :: blat_in
     integer :: i, j, nlon_in, nlat_in
     real(r4_kind)    :: tpi
     integer, parameter    :: kindl = r4_kind !< real size at compile time
!-----------------------------------------------------------------------
   call horiz_interp_init

   tpi = 2.0_kindl * real(pi, r4_kind)
   nlon_in = size(data_in,1)
   nlat_in = size(data_in,2)

   do i = 1, nlon_in+1
      blon_in(i) = wb + real(i-1, r4_kind)*dx
   enddo
      if (abs(blon_in(nlon_in+1)-blon_in(1)-tpi) < epsilon(blon_in)) &
              blon_in(nlon_in+1)=blon_in(1)+tpi

   do j = 2, nlat_in
      blat_in(j) = sb + real(j-1, r4_kind)*dy
   enddo
      blat_in(1)         = -0.5_kindl * real(pi, r4_kind)
      blat_in(nlat_in+1) =  0.5_kindl * real(pi, r4_kind)


   call horiz_interp_solo_1d (data_in, blon_in, blat_in,    &
                              lon_out, lat_out, data_out,   &
                              verbose, mask_in, mask_out    )

!-----------------------------------------------------------------------

 end subroutine horiz_interp_solo_old_r4

!#######################################################################


 !####################################################################
 function is_lat_lon_r4(lon, lat)
    real(r4_kind), dimension(:,:), intent(in) :: lon, lat
    logical                          :: is_lat_lon_r4
    integer                          :: i, j, nlon, nlat, num

    is_lat_lon_r4 = .true.
    nlon = size(lon,1)
    nlat = size(lon,2)
    LOOP_LAT: do j = 1, nlat
       do i = 2, nlon
          if(lat(i,j) .NE. lat(1,j)) then
             is_lat_lon_r4 = .false.
             exit LOOP_LAT
          end if
       end do
    end do LOOP_LAT

    if(is_lat_lon_r4) then
       LOOP_LON: do i = 1, nlon
          do j = 2, nlat
             if(lon(i,j) .NE. lon(i,1)) then
                is_lat_lon_r4 = .false.
                exit LOOP_LON
             end if
          end do
       end do LOOP_LON
    end if

    num = 0
    if(is_lat_lon_r4) num = 1
    call mpp_min(num)
    if(num == 1) then
       is_lat_lon_r4 = .true.
    else
       is_lat_lon_r4 = .false.
    end if

    return
 end function is_lat_lon_r4

 !> Subroutine for reading a weight file and use it to fill in the horiz interp type
!! for the bilinear interpolation method.
 subroutine horiz_interp_read_weights_r4(Interp, weight_filename, lon_out, lat_out, lon_in, lat_in, &
                                       weight_file_source, interp_method, isw, iew, jsw, jew, nglon, nglat)
   type(horiz_interp_type), intent(inout) :: Interp             !< Horiz interp time to fill
   character(len=*),        intent(in)    :: weight_filename    !< Name of the weight file
   real(r4_kind),      intent(in)    :: lat_out(:,:)       !< Output (model) latitude
   real(r4_kind),      intent(in)    :: lon_out(:,:)       !< Output (model) longitude
   real(r4_kind),      intent(in)    :: lat_in(:)          !< Input (data) latitude
   real(r4_kind),      intent(in)    :: lon_in(:)          !< Input (data) longitude
   character(len=*),        intent(in)    :: weight_file_source !< Source of the weight file
   character(len=*),        intent(in)    :: interp_method      !< The interp method to use
   integer,                 intent(in)    :: isw, iew, jsw, jew !< Starting and ending indices of the compute domain
   integer,                 intent(in)    :: nglon              !< Number of longitudes in the global domain
   integer,                 intent(in)    :: nglat              !< Number of latitudes in the globl domain

   integer                           :: i, j                 !< For do loops
   integer                           :: nlon_in              !< Number of longitude in the data
   integer                           :: nlat_in              !< Number of latitude in the data grid
   real(r4_kind), allocatable   :: lon_src_1d(:)        !< Center points of the longitude data grid
   real(r4_kind), allocatable   :: lat_src_1d(:)        !< Center points of the lattiude data grid
   integer, parameter                :: kindl = r4_kind !< real kind size currently compiling

   select case (trim(interp_method))
   case ("bilinear")
     !! This is to reproduce the behavior in horiz_interp_new
     !! The subroutine assumes that the data grid (lon_in, lat_in) are
     !! the edges and not the centers.
     !! Data_override passes in the edges, which are calculated using the axis_edges subroutine
     nlon_in  = size(lon_in(:))-1;  nlat_in  = size(lat_in(:))-1
     allocate(lon_src_1d(nlon_in), lat_src_1d(nlat_in))
     do i = 1, nlon_in
       lon_src_1d(i) = (lon_in(i) + lon_in(i+1)) * 0.5_kindl
     enddo
     do j = 1, nlat_in
       lat_src_1d(j) = (lat_in(j) + lat_in(j+1)) * 0.5_kindl
     enddo

     call horiz_interp_read_weights_bilinear(Interp, weight_filename, lon_out, lat_out, &
      lon_src_1d, lat_src_1d, weight_file_source, interp_method, &
      isw, iew, jsw, jew, nglon, nglat)
     deallocate(lon_src_1d,lat_src_1d)
   case default
     call mpp_error(FATAL, "Reading weight from file is not supported for the "//&
      trim(interp_method)//" method. It is currently only supported for bilinear")
   end select
 end subroutine horiz_interp_read_weights_r4
# 64 "horiz_interp/include/horiz_interp_r4.fh" 2

# 318 "horiz_interp/horiz_interp.F90" 2

# 1 "horiz_interp/include/horiz_interp_r8.fh" 1
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














































# 1 "horiz_interp/include/horiz_interp.inc" 1
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
  !> @brief Creates a 1D @ref horiz_interp_type with the given parameters
  subroutine horiz_interp_new_1d_r8 (Interp, lon_in, lat_in, lon_out, lat_out, verbose, &
                                  interp_method, num_nbrs, max_dist, src_modulo,     &
                                  grid_at_center, mask_in, mask_out)

    !-----------------------------------------------------------------------
    type(horiz_interp_type), intent(inout)        :: Interp
    real(r8_kind), intent(in),  dimension(:)               :: lon_in , lat_in
    real(r8_kind), intent(in),  dimension(:)               :: lon_out, lat_out
    integer, intent(in),                 optional :: verbose
    character(len=*), intent(in),        optional :: interp_method
    integer, intent(in),                 optional :: num_nbrs
    real(r8_kind),    intent(in),                 optional :: max_dist
    logical, intent(in),                 optional :: src_modulo
    logical, intent(in),                 optional :: grid_at_center
    real(r8_kind), intent(in), dimension(:,:),    optional :: mask_in  !< dummy variable
    real(r8_kind), intent(inout),dimension(:,:),  optional :: mask_out !< dummy variable
    !-----------------------------------------------------------------------
    real(r8_kind), dimension(:,:), allocatable :: lon_src, lat_src, lon_dst, lat_dst
    real(r8_kind), dimension(:),   allocatable :: lon_src_1d, lat_src_1d, lon_dst_1d, lat_dst_1d
    integer                           :: i, j, nlon_in, nlat_in, nlon_out, nlat_out
    logical                           :: center
    character(len=40)                 :: method
    integer, parameter                :: kindl = r8_kind !> real kind size currently compiling
    !-----------------------------------------------------------------------
    call horiz_interp_init

    method = 'conservative'
    if(present(interp_method)) method = interp_method

    select case (trim(method))
    case ("conservative")
       Interp%interp_method = CONSERVE
       call horiz_interp_conserve_new ( Interp, lon_in, lat_in, lon_out, lat_out, verbose)
    case ("bilinear")
       Interp%interp_method = BILINEAR
       center = .false.
       if(present(grid_at_center) ) center = grid_at_center
       if(center) then
          nlon_out = size(lon_out(:)); nlat_out = size(lat_out(:))
          allocate(lon_dst(nlon_out,nlat_out), lat_dst(nlon_out,nlat_out))
          do i = 1, nlon_out
             lon_dst(i,:) = lon_out(i)
          enddo
          do j = 1, nlat_out
             lat_dst(:,j) = lat_out(j)
          enddo

          call horiz_interp_bilinear_new ( Interp, lon_in, lat_in, lon_dst, lat_dst, &
               verbose, src_modulo)
          deallocate(lon_dst, lat_dst)
       else
          nlon_in  = size(lon_in(:))-1;  nlat_in  = size(lat_in(:))-1
          nlon_out = size(lon_out(:))-1; nlat_out = size(lat_out(:))-1
          allocate(lon_src_1d(nlon_in), lat_src_1d(nlat_in))
          allocate(lon_dst(nlon_out,nlat_out), lat_dst(nlon_out,nlat_out))
          do i = 1, nlon_in
             lon_src_1d(i) = (lon_in(i) + lon_in(i+1)) * 0.5_kindl
          enddo
          do j = 1, nlat_in
             lat_src_1d(j) = (lat_in(j) + lat_in(j+1)) * 0.5_kindl
          enddo
          do i = 1, nlon_out
             lon_dst(i,:) = (lon_out(i) + lon_out(i+1)) * 0.5_kindl
          enddo
          do j = 1, nlat_out
             lat_dst(:,j) = (lat_out(j) + lat_out(j+1)) * 0.5_kindl
          enddo
          call horiz_interp_bilinear_new ( Interp, lon_src_1d, lat_src_1d, lon_dst, lat_dst, &
               verbose, src_modulo)
          deallocate(lon_src_1d, lat_src_1d, lon_dst, lat_dst)
       endif
    case ("bicubic")
       Interp%interp_method = BICUBIC
       center = .false.
       if(present(grid_at_center) ) center = grid_at_center
       !No need to expand to 2d, horiz_interp_bicubic_new does 1d-1d
       if(center) then
          call horiz_interp_bicubic_new ( Interp, lon_in, lat_in, lon_out, lat_out, &
            verbose, src_modulo)
       else
          nlon_in  = size(lon_in(:))-1;  nlat_in  = size(lat_in(:))-1
          nlon_out = size(lon_out(:))-1; nlat_out = size(lat_out(:))-1
          allocate(lon_src_1d(nlon_in), lat_src_1d(nlat_in))
          allocate(lon_dst_1d(nlon_out), lat_dst_1d(nlat_out))
          do i = 1, nlon_in
             lon_src_1d(i) = (lon_in(i) + lon_in(i+1)) * 0.5_kindl
          enddo
          do j = 1, nlat_in
             lat_src_1d(j) = (lat_in(j) + lat_in(j+1)) * 0.5_kindl
          enddo
          do i = 1, nlon_out
             lon_dst_1d(i) = (lon_out(i) + lon_out(i+1)) * 0.5_kindl
          enddo
          do j = 1, nlat_out
             lat_dst_1d(j) = (lat_out(j) + lat_out(j+1)) * 0.5_kindl
          enddo
          call horiz_interp_bicubic_new ( Interp, lon_src_1d, lat_src_1d, lon_dst_1d, lat_dst_1d, &
               verbose, src_modulo)
          deallocate(lon_src_1d, lat_src_1d, lon_dst_1d, lat_dst_1d)
       endif
    case ("spherical")
       Interp%interp_method = SPHERICAL
       nlon_in  = size(lon_in(:));   nlat_in  = size(lat_in(:))
       nlon_out  = size(lon_out(:)); nlat_out = size(lat_out(:))
       allocate(lon_src(nlon_in,nlat_in), lat_src(nlon_in,nlat_in))
       allocate(lon_dst(nlon_out,nlat_out), lat_dst(nlon_out,nlat_out))
       do i = 1, nlon_in
          lon_src(i,:) = lon_in(i)
       enddo
       do j = 1, nlat_in
          lat_src(:,j) = lat_in(j)
       enddo
       do i = 1, nlon_out
          lon_dst(i,:) = lon_out(i)
       enddo
       do j = 1, nlat_out
          lat_dst(:,j) = lat_out(j)
       enddo
       call horiz_interp_spherical_new ( Interp, lon_src, lat_src, lon_dst, lat_dst, &
            num_nbrs, max_dist, src_modulo)
       deallocate(lon_src, lat_src, lon_dst, lat_dst)
    case default
       call mpp_error(FATAL,'horiz_interp_mod: interp_method should be conservative, bilinear, bicubic, spherical')
    end select

    !-----------------------------------------------------------------------
    Interp% horizInterpReals8_type % is_allocated = .true.
    Interp%I_am_initialized = .true.

  end subroutine horiz_interp_new_1d_r8

!#######################################################################

 subroutine horiz_interp_new_1d_src_r8 (Interp, lon_in, lat_in, lon_out, lat_out,   &
                                     verbose, interp_method, num_nbrs, max_dist, &
                                     src_modulo, grid_at_center, mask_in, mask_out, is_latlon_out )

   type(horiz_interp_type), intent(inout)        :: Interp
   real(r8_kind), intent(in),  dimension(:)               :: lon_in , lat_in
   real(r8_kind), intent(in),  dimension(:,:)             :: lon_out, lat_out
   integer, intent(in),                 optional :: verbose
   character(len=*), intent(in),        optional :: interp_method
   integer, intent(in),                 optional :: num_nbrs  !< minimum number of neighbors
   real(r8_kind),    intent(in),                 optional :: max_dist
   logical, intent(in),                 optional :: src_modulo
   logical, intent(in),                 optional :: grid_at_center
   real(r8_kind), intent(in), dimension(:,:),    optional :: mask_in
   real(r8_kind), intent(out),dimension(:,:),    optional :: mask_out
   logical, intent(in),                 optional :: is_latlon_out

   real(r8_kind), dimension(:,:), allocatable :: lon_src, lat_src
   real(r8_kind), dimension(:),   allocatable :: lon_src_1d, lat_src_1d
   integer                           :: i, j, nlon_in, nlat_in
   character(len=40)                 :: method
   logical                           :: center
   logical                           :: dst_is_latlon
   integer, parameter                :: kindl = r8_kind !< real kind size currently compiling
   !-----------------------------------------------------------------------
   call horiz_interp_init

   method = 'conservative'
   if(present(interp_method)) method = interp_method

   select case (trim(method))
   case ("conservative")
      Interp%interp_method = CONSERVE
      !--- check to see if the source grid is regular lat-lon grid or not.
      if(PRESENT(is_latlon_out)) then
         dst_is_latlon = is_latlon_out
      else
         dst_is_latlon = is_lat_lon(lon_out, lat_out)
      end if
      if(dst_is_latlon ) then
         if(present(mask_in)) then
            if ( ANY(mask_in < -.0001_kindl) .or. ANY(mask_in > 1.0001_kindl)) &
               call mpp_error(FATAL, &
                  'horiz_interp_conserve_new_1d_src(horiz_interp_conserve_mod): input mask not between 0,1')
            allocate(Interp%horizInterpReals8_type%mask_in(size(mask_in,1), size(mask_in,2)) )
            Interp%horizInterpReals8_type%mask_in = mask_in
         end if
         call horiz_interp_conserve_new ( Interp, lon_in, lat_in, lon_out(:,1), lat_out(1,:), &
              verbose=verbose )
      else
         call horiz_interp_conserve_new ( Interp, lon_in, lat_in, lon_out, lat_out, &
              verbose=verbose, mask_in=mask_in, mask_out=mask_out )
      end if
   case ("bilinear")
      Interp%interp_method = BILINEAR
      center = .false.
      if(present(grid_at_center) ) center = grid_at_center
      if(center) then
         call horiz_interp_bilinear_new ( Interp, lon_in, lat_in, lon_out, lat_out, &
              verbose, src_modulo )
      else
         nlon_in  = size(lon_in(:))-1;  nlat_in  = size(lat_in(:))-1
         allocate(lon_src_1d(nlon_in), lat_src_1d(nlat_in))
         do i = 1, nlon_in
            lon_src_1d(i) = (lon_in(i) + lon_in(i+1)) * 0.5_kindl
         enddo
         do j = 1, nlat_in
            lat_src_1d(j) = (lat_in(j) + lat_in(j+1)) * 0.5_kindl
         enddo
         call horiz_interp_bilinear_new ( Interp, lon_src_1d, lat_src_1d, lon_out, lat_out, &
              verbose, src_modulo )
         deallocate(lon_src_1d,lat_src_1d)
      endif
   case ("bicubic")
      Interp%interp_method = BICUBIC
      center = .false.
      if(present(grid_at_center) ) center = grid_at_center
      if(center) then
        call horiz_interp_bicubic_new ( Interp, lon_in, lat_in, lon_out, lat_out, &
              verbose, src_modulo )
      else
         nlon_in  = size(lon_in(:))-1;  nlat_in  = size(lat_in(:))-1
         allocate(lon_src_1d(nlon_in), lat_src_1d(nlat_in))
         do i = 1, nlon_in
            lon_src_1d(i) = (lon_in(i) + lon_in(i+1)) * 0.5_kindl
         enddo
         do j = 1, nlat_in
            lat_src_1d(j) = (lat_in(j) + lat_in(j+1)) * 0.5_kindl
         enddo
           call horiz_interp_bicubic_new ( Interp, lon_src_1d, lat_src_1d, lon_out, lat_out, &
              verbose, src_modulo )
         deallocate(lon_src_1d,lat_src_1d)
      endif
   case ("spherical")
      Interp%interp_method = SPHERICAL
      nlon_in  = size(lon_in(:));  nlat_in  = size(lat_in(:))
      allocate(lon_src(nlon_in,nlat_in), lat_src(nlon_in,nlat_in))
      do i = 1, nlon_in
         lon_src(i,:) = lon_in(i)
      enddo
      do j = 1, nlat_in
         lat_src(:,j) = lat_in(j)
      enddo
      call horiz_interp_spherical_new ( Interp, lon_src, lat_src, lon_out, lat_out, &
           num_nbrs, max_dist, src_modulo)
      deallocate(lon_src, lat_src)
   case default
      call mpp_error(FATAL,'interp_method should be conservative, bilinear, bicubic, spherical')
   end select

   !-----------------------------------------------------------------------
   Interp% horizInterpReals8_type % is_allocated = .true.
   Interp%I_am_initialized = .true.

 end subroutine horiz_interp_new_1d_src_r8

!#######################################################################

 subroutine horiz_interp_new_2d_r8 (Interp, lon_in, lat_in, lon_out, lat_out,   &
                                 verbose, interp_method, num_nbrs, max_dist, &
                                 src_modulo, mask_in, mask_out, is_latlon_in, is_latlon_out  )
 type(horiz_interp_type), intent(inout)     :: Interp
 real(r8_kind), intent(in),  dimension(:,:)          :: lon_in , lat_in
 real(r8_kind), intent(in),  dimension(:,:)          :: lon_out, lat_out
 integer, intent(in),              optional :: verbose
 character(len=*), intent(in),     optional :: interp_method
 integer, intent(in),              optional :: num_nbrs
 real(r8_kind),    intent(in),              optional :: max_dist
 logical, intent(in),              optional :: src_modulo
 real(r8_kind), intent(in), dimension(:,:), optional :: mask_in
 real(r8_kind), intent(out),dimension(:,:), optional :: mask_out
 logical, intent(in),              optional :: is_latlon_in, is_latlon_out
 logical           :: src_is_latlon, dst_is_latlon
 character(len=40) :: method
 integer, parameter                :: kindl = r8_kind !< real kind size currently compiling
!-----------------------------------------------------------------------
   call horiz_interp_init

   method = 'bilinear'
   if(present(interp_method)) method = interp_method

   select case (trim(method))
   case ("conservative")
      Interp%interp_method = CONSERVE
      if(PRESENT(is_latlon_in)) then
         src_is_latlon = is_latlon_in
      else
         src_is_latlon = is_lat_lon(lon_in, lat_in)
      end if
      if(PRESENT(is_latlon_out)) then
         dst_is_latlon = is_latlon_out
      else
         dst_is_latlon = is_lat_lon(lon_out, lat_out)
      end if
      if(src_is_latlon .AND. dst_is_latlon) then
         if(present(mask_in)) then
            if ( ANY(mask_in < -0.0001_kindl) .or. ANY(mask_in > 1.0001_kindl)) then
                call mpp_error(FATAL, 'horiz_interp_conserve_new_2d(horiz_interp_conserve_mod):' // &
                                       ' input mask not between 0,1')
            endif
            allocate(Interp%horizInterpReals8_type%mask_in(size(mask_in,1), size(mask_in,2)) )
            Interp%horizInterpReals8_type%mask_in = mask_in
         end if
         call horiz_interp_conserve_new ( Interp, lon_in(:,1), lat_in(1,:), lon_out(:,1), lat_out(1,:), &
              verbose=verbose )
      else if(src_is_latlon) then
         call horiz_interp_conserve_new ( Interp, lon_in(:,1), lat_in(1,:), lon_out, lat_out, &
              verbose=verbose, mask_in=mask_in, mask_out=mask_out )
      else if(dst_is_latlon) then
         call horiz_interp_conserve_new ( Interp, lon_in, lat_in, lon_out(:,1), lat_out(1,:), &
              verbose=verbose, mask_in=mask_in, mask_out=mask_out )
      else
         call horiz_interp_conserve_new ( Interp, lon_in, lat_in, lon_out, lat_out, &
              verbose=verbose, mask_in=mask_in, mask_out=mask_out )
      end if

   case ("spherical")
      Interp%interp_method = SPHERICAL
      call horiz_interp_spherical_new ( Interp, lon_in, lat_in, lon_out, lat_out, &
                                    num_nbrs, max_dist, src_modulo )
   case ("bilinear")
      Interp%interp_method = BILINEAR
      call horiz_interp_bilinear_new ( Interp, lon_in, lat_in, lon_out, lat_out, &
                                        verbose, src_modulo )
   case default
      call mpp_error(FATAL,'when source grid are 2d, interp_method should be spherical or bilinear')
   end select

   Interp% horizInterpReals8_type % is_allocated = .true.
   Interp%I_am_initialized = .true.

 end subroutine horiz_interp_new_2d_r8

!#######################################################################
 subroutine horiz_interp_new_1d_dst_r8 (Interp, lon_in, lat_in, lon_out, lat_out,   &
      verbose, interp_method, num_nbrs, max_dist, src_modulo, mask_in, mask_out, is_latlon_in )
   type(horiz_interp_type), intent(inout)     :: Interp
   real(r8_kind), intent(in),  dimension(:,:)          :: lon_in , lat_in
   real(r8_kind), intent(in),  dimension(:)            :: lon_out, lat_out
   integer, intent(in),              optional :: verbose
   character(len=*), intent(in),     optional :: interp_method
   integer, intent(in),              optional :: num_nbrs
   real(r8_kind),    intent(in),              optional :: max_dist
   logical, intent(in),              optional :: src_modulo
   real(r8_kind), intent(in), dimension(:,:), optional :: mask_in
   real(r8_kind), intent(out),dimension(:,:), optional :: mask_out
   logical, intent(in),              optional :: is_latlon_in

   character(len=40) :: method
   integer, parameter                :: kindl = r8_kind !< real kind size currently compiling
   !-------------some local variables-----------------------------------------------
   integer                           :: i, j, nlon_out, nlat_out
   real(r8_kind), dimension(:,:), allocatable :: lon_dst, lat_dst
   logical                           :: src_is_latlon
   !-----------------------------------------------------------------------
   call horiz_interp_init

   method = 'bilinear'
   if(present(interp_method)) method = interp_method

   nlon_out = size(lon_out(:)); nlat_out = size(lat_out(:))
   allocate(lon_dst(nlon_out,nlat_out), lat_dst(nlon_out,nlat_out))
   do i = 1, nlon_out
      lon_dst(i,:) = lon_out(i)
   enddo
   do j = 1, nlat_out
      lat_dst(:,j) = lat_out(j)
   enddo

   select case (trim(method))
   case ("conservative")
      Interp%interp_method = CONSERVE
      if(PRESENT(is_latlon_in)) then
         src_is_latlon = is_latlon_in
      else
         src_is_latlon = is_lat_lon(lon_in, lat_in)
      end if

      if(src_is_latlon) then
         if(present(mask_in)) then
            if ( ANY(mask_in < -0.0001_kindl) .or. ANY(mask_in > 1.0001_kindl)) &
                 call mpp_error(FATAL, &
                     'horiz_interp_conserve_new_1d_dst(horiz_interp_conserve_mod): input mask not between 0,1')
            allocate(Interp%horizInterpReals8_type%mask_in(size(mask_in,1), size(mask_in,2)) )
            Interp%horizInterpReals8_type%mask_in = mask_in
         end if
         call horiz_interp_conserve_new ( Interp, lon_in(:,1), lat_in(1,:), lon_out, lat_out, &
              verbose=verbose)
      else
         call horiz_interp_conserve_new ( Interp, lon_in, lat_in, lon_out, lat_out, &
              verbose=verbose, mask_in=mask_in, mask_out=mask_out )
      end if
   case ("bilinear")
      Interp%interp_method = BILINEAR
      call horiz_interp_bilinear_new ( Interp, lon_in, lat_in, lon_dst, lat_dst, &
           verbose, src_modulo )
   case ("spherical")
      Interp%interp_method = SPHERICAL
      call horiz_interp_spherical_new ( Interp, lon_in, lat_in, lon_dst, lat_dst, &
           num_nbrs, max_dist, src_modulo)
   case default
      call mpp_error(FATAL,'when source grid are 2d, interp_method should be spherical or bilinear')
   end select

   deallocate(lon_dst,lat_dst)

   Interp% horizInterpReals8_type % is_allocated = .true.
   Interp%I_am_initialized = .true.

 end subroutine horiz_interp_new_1d_dst_r8

!#######################################################################

 subroutine horiz_interp_base_2d_r8 ( Interp, data_in, data_out, verbose, &
                                   mask_in, mask_out, missing_value, missing_permit, &
                                   err_msg, new_missing_handle )
!-----------------------------------------------------------------------
   type (horiz_interp_type), intent(in) :: Interp
      real(r8_kind), intent(in),  dimension(:,:) :: data_in
      real(r8_kind), intent(out), dimension(:,:) :: data_out
   integer, intent(in),                   optional :: verbose
      real(r8_kind), intent(in),   dimension(:,:), optional :: mask_in
      real(r8_kind), intent(out),  dimension(:,:), optional :: mask_out
      real(r8_kind), intent(in),                   optional :: missing_value
      integer, intent(in),                optional :: missing_permit
   character(len=*), intent(out),         optional :: err_msg
      logical, intent(in),                optional :: new_missing_handle
!-----------------------------------------------------------------------
   if(present(err_msg)) err_msg = ''
   if(.not.Interp%I_am_initialized) then
     if(fms_error_handler('horiz_interp','The horiz_interp_type variable is not initialized',err_msg)) return
   endif

   select case(Interp%interp_method)
   case(CONSERVE)
      call horiz_interp_conserve(Interp,data_in, data_out, verbose, mask_in, mask_out)
   case(BILINEAR)
      call horiz_interp_bilinear(Interp,data_in, data_out, verbose, mask_in, mask_out, &
                             missing_value, missing_permit, new_missing_handle )
   case(BICUBIC)
      call horiz_interp_bicubic(Interp,data_in, data_out, verbose, mask_in, mask_out, &
                             missing_value, missing_permit )
   case(SPHERICAL)
      call horiz_interp_spherical(Interp,data_in, data_out, verbose, mask_in, mask_out, &
                             missing_value )
   case default
      call mpp_error(FATAL,'interp_method should be conservative, bilinear, bicubic, spherical')
   end select

   return

 end subroutine horiz_interp_base_2d_r8

!#######################################################################

 !> Overload of interface horiz_interp_base_2d_r8
 !! uses 3d arrays for data and mask
 !! this allows for multiple interpolations with one call
 subroutine horiz_interp_base_3d_r8 ( Interp, data_in, data_out, verbose, mask_in, mask_out, &
      missing_value, missing_permit, err_msg  )
   !-----------------------------------------------------------------------
   !   overload of interface horiz_interp_base_2d_r8
   !   uses 3d arrays for data and mask
   !   this allows for multiple interpolations with one call
   !-----------------------------------------------------------------------
   type (horiz_interp_type), intent(in)           :: Interp
   real(r8_kind), intent(in),  dimension(:,:,:)            :: data_in
   real(r8_kind), intent(out), dimension(:,:,:)            :: data_out
   integer, intent(in),                  optional :: verbose
   real(r8_kind), intent(in),   dimension(:,:,:), optional :: mask_in
   real(r8_kind), intent(out),  dimension(:,:,:), optional :: mask_out
   real(r8_kind), intent(in),                     optional :: missing_value
   integer, intent(in),                  optional :: missing_permit
   character(len=*), intent(out),        optional :: err_msg
   !-----------------------------------------------------------------------
   integer :: n

   if(present(err_msg)) err_msg = ''
   if(.not.Interp%I_am_initialized) then
     if(fms_error_handler('horiz_interp','The horiz_interp_type variable is not initialized',err_msg)) return
   endif

   do n = 1, size(data_in,3)
      if (present(mask_in))then
         if(present(mask_out)) then
            call horiz_interp( Interp, data_in(:,:,n), data_out(:,:,n), &
                 verbose, mask_in(:,:,n), mask_out(:,:,n), &
                 missing_value, missing_permit )
         else
            call horiz_interp( Interp, data_in(:,:,n), data_out(:,:,n), &
                 verbose, mask_in(:,:,n), missing_value = missing_value,  &
                 missing_permit = missing_permit )
         endif
      else
         if(present(mask_out)) then
            call horiz_interp( Interp, data_in(:,:,n), data_out(:,:,n), &
                 verbose, mask_out=mask_out(:,:,n), missing_value = missing_value,  &
                 missing_permit = missing_permit )
         else
            call horiz_interp( Interp, data_in(:,:,n), data_out(:,:,n), &
                 verbose, missing_value = missing_value,  &
                 missing_permit = missing_permit )
         endif
     endif
   enddo

   return
!-----------------------------------------------------------------------
 end subroutine horiz_interp_base_3d_r8

!#######################################################################

!> Interpolates from a rectangular grid to rectangular grid.
!! interp_method can be the value conservative, bilinear or spherical.
!! horiz_interp_new don't need to be called before calling this routine.
 subroutine horiz_interp_solo_1d_r8 ( data_in, lon_in, lat_in, lon_out, lat_out,    &
                                   data_out, verbose, mask_in, mask_out,         &
                                   interp_method, missing_value, missing_permit, &
                                   num_nbrs, max_dist,src_modulo, grid_at_center  )
!-----------------------------------------------------------------------
      real(r8_kind), intent(in),  dimension(:,:) :: data_in
      real(r8_kind), intent(in),  dimension(:)   :: lon_in , lat_in
      real(r8_kind), intent(in),  dimension(:)   :: lon_out, lat_out
      real(r8_kind), intent(out), dimension(:,:) :: data_out
   integer, intent(in),                   optional :: verbose
      real(r8_kind), intent(in),   dimension(:,:), optional :: mask_in
      real(r8_kind), intent(out),  dimension(:,:), optional :: mask_out
   character(len=*), intent(in),          optional :: interp_method
      real(r8_kind), intent(in),                   optional :: missing_value
   integer, intent(in),                   optional :: missing_permit
   integer, intent(in),                   optional :: num_nbrs
      real(r8_kind), intent(in),                   optional :: max_dist
   logical, intent(in),                   optional :: src_modulo
   logical, intent(in),                   optional :: grid_at_center
!-----------------------------------------------------------------------
    type (horiz_interp_type) :: Interp
!-----------------------------------------------------------------------
    call horiz_interp_init

    call horiz_interp_new ( Interp, lon_in, lat_in, lon_out, lat_out, verbose, &
                             interp_method, num_nbrs, max_dist, src_modulo, grid_at_center )

    call horiz_interp ( Interp, data_in, data_out, verbose,   &
                        mask_in, mask_out, missing_value, missing_permit )

    call horiz_interp_del ( Interp )
!-----------------------------------------------------------------------

 end subroutine horiz_interp_solo_1d_r8

!#######################################################################

!> Interpolates from a uniformly spaced grid to any output grid.
!! interp_method can be the value "onservative","bilinear" or "spherical".
!! horiz_interp_new don't need to be called before calling this routine.
 subroutine horiz_interp_solo_1d_src_r8 ( data_in, lon_in, lat_in, lon_out, lat_out,    &
                                       data_out, verbose, mask_in, mask_out,         &
                                       interp_method, missing_value, missing_permit, &
                                       num_nbrs, max_dist, src_modulo, grid_at_center )
!-----------------------------------------------------------------------
      real(r8_kind), intent(in),  dimension(:,:) :: data_in
      real(r8_kind), intent(in),  dimension(:)   :: lon_in , lat_in
      real(r8_kind), intent(in),  dimension(:,:) :: lon_out, lat_out
      real(r8_kind), intent(out), dimension(:,:) :: data_out
   integer, intent(in),                   optional :: verbose
      real(r8_kind), intent(in),   dimension(:,:), optional :: mask_in
      real(r8_kind), intent(out),  dimension(:,:), optional :: mask_out
   character(len=*), intent(in),          optional :: interp_method
      real(r8_kind), intent(in),                   optional :: missing_value
   integer, intent(in),                   optional :: missing_permit
   integer, intent(in),                   optional :: num_nbrs
      real(r8_kind), intent(in),                   optional :: max_dist
   logical, intent(in),                   optional :: src_modulo
   logical, intent(in),                   optional :: grid_at_center

!-----------------------------------------------------------------------
   type (horiz_interp_type) :: Interp
   logical                  :: dst_is_latlon
   character(len=128)       :: method
!-----------------------------------------------------------------------
    call horiz_interp_init
    method = 'conservative'
    if(present(interp_method)) method = interp_method
    dst_is_latlon = .true.
    if(trim(method) == 'conservative') dst_is_latlon = is_lat_lon(lon_out, lat_out)

    if(dst_is_latlon) then
       call horiz_interp_new ( Interp, lon_in, lat_in, lon_out, lat_out, verbose, &
                               interp_method, num_nbrs, max_dist, src_modulo,    &
                               grid_at_center, is_latlon_out = dst_is_latlon )
       call horiz_interp ( Interp, data_in, data_out, verbose,   &
                           mask_in, mask_out, missing_value, missing_permit )
    else
       call horiz_interp_new ( Interp, lon_in, lat_in, lon_out, lat_out, verbose, &
                               interp_method, num_nbrs, max_dist, src_modulo,    &
                               grid_at_center, mask_in, mask_out, is_latlon_out = dst_is_latlon)

       call horiz_interp ( Interp, data_in, data_out, verbose,   &
                           missing_value=missing_value, missing_permit=missing_permit )
    end if

    call horiz_interp_del ( Interp )

!-----------------------------------------------------------------------

 end subroutine horiz_interp_solo_1d_src_r8


!#######################################################################

!> Interpolates from any grid to any grid. interp_method should be "spherical"
!! horiz_interp_new don't need to be called before calling this routine.
 subroutine horiz_interp_solo_2d_r8 ( data_in, lon_in, lat_in, lon_out, lat_out, data_out, &
                                   verbose, mask_in, mask_out, interp_method, missing_value,&
                                   missing_permit, num_nbrs, max_dist, src_modulo  )
!-----------------------------------------------------------------------
      real(r8_kind), intent(in),  dimension(:,:) :: data_in
      real(r8_kind), intent(in),  dimension(:,:) :: lon_in , lat_in
      real(r8_kind), intent(in),  dimension(:,:) :: lon_out, lat_out
      real(r8_kind), intent(out), dimension(:,:) :: data_out
   integer, intent(in),                   optional :: verbose
      real(r8_kind), intent(in),   dimension(:,:), optional :: mask_in
      real(r8_kind), intent(out),  dimension(:,:), optional :: mask_out
   character(len=*), intent(in),          optional :: interp_method
      real(r8_kind), intent(in),                   optional :: missing_value
   integer, intent(in),                   optional :: missing_permit
   integer, intent(in),                   optional :: num_nbrs
      real(r8_kind), intent(in),                   optional :: max_dist
   logical, intent(in),                   optional :: src_modulo
!-----------------------------------------------------------------------
   type (horiz_interp_type) :: Interp
   logical                  :: dst_is_latlon, src_is_latlon
   character(len=128)       :: method
!-----------------------------------------------------------------------
    call horiz_interp_init

    method = 'conservative'
    if(present(interp_method)) method = interp_method
    dst_is_latlon = .true.
    src_is_latlon = .true.
    if(trim(method) == 'conservative') then
       dst_is_latlon = is_lat_lon(lon_out, lat_out)
       src_is_latlon = is_lat_lon(lon_in, lat_in)
    end if

    if(dst_is_latlon .and. src_is_latlon) then
       call horiz_interp_new ( Interp, lon_in, lat_in, lon_out, lat_out, verbose, &
                               interp_method, num_nbrs, max_dist, src_modulo,    &
                               is_latlon_in=dst_is_latlon, is_latlon_out = dst_is_latlon )
       call horiz_interp ( Interp, data_in, data_out, verbose,   &
                           mask_in, mask_out, missing_value, missing_permit )
    else
       call horiz_interp_new ( Interp, lon_in, lat_in, lon_out, lat_out, verbose, &
                               interp_method, num_nbrs, max_dist, src_modulo,    &
                               mask_in, mask_out, &
                               is_latlon_in=dst_is_latlon, is_latlon_out = dst_is_latlon)
       call horiz_interp ( Interp, data_in, data_out, verbose,   &
                           missing_value=missing_value, missing_permit=missing_permit )
    end if

    call horiz_interp_del ( Interp )

!-----------------------------------------------------------------------

 end subroutine horiz_interp_solo_2d_r8

!#######################################################################

!>   interpolates from any grid to rectangular longitude/latitude grid.
!!   interp_method should be "spherical".
!!   horiz_interp_new don't need to be called before calling this routine.
 subroutine horiz_interp_solo_1d_dst_r8 ( data_in, lon_in, lat_in, lon_out, lat_out, data_out,    &
                                       verbose, mask_in, mask_out,interp_method,missing_value, &
                                       missing_permit,  num_nbrs, max_dist, src_modulo)
!-----------------------------------------------------------------------
      real(r8_kind), intent(in),  dimension(:,:) :: data_in
      real(r8_kind), intent(in),  dimension(:,:) :: lon_in , lat_in
      real(r8_kind), intent(in),  dimension(:)   :: lon_out, lat_out
      real(r8_kind), intent(out), dimension(:,:) :: data_out
   integer, intent(in),                   optional :: verbose
      real(r8_kind), intent(in),   dimension(:,:), optional :: mask_in
      real(r8_kind), intent(out),  dimension(:,:), optional :: mask_out
   character(len=*), intent(in),          optional :: interp_method
      real(r8_kind), intent(in),                   optional :: missing_value
   integer, intent(in),                   optional :: missing_permit
   integer, intent(in),                   optional :: num_nbrs
      real(r8_kind), intent(in),                   optional :: max_dist
   logical, intent(in),                   optional :: src_modulo
!-----------------------------------------------------------------------
   type (horiz_interp_type) :: Interp
   logical                  :: src_is_latlon
   character(len=128)       :: method
!-----------------------------------------------------------------------
    call horiz_interp_init

    method = 'conservative'
    if(present(interp_method)) method = interp_method
    src_is_latlon = .true.
    if(trim(method) == 'conservative') src_is_latlon = is_lat_lon(lon_in, lat_in)

    if(src_is_latlon) then
       call horiz_interp_new ( Interp, lon_in, lat_in, lon_out, lat_out, verbose, &
                               interp_method, num_nbrs, max_dist, src_modulo,    &
                               is_latlon_in = src_is_latlon )
       call horiz_interp ( Interp, data_in, data_out, verbose,   &
                           mask_in, mask_out, missing_value, missing_permit )
    else
       call horiz_interp_new ( Interp, lon_in, lat_in, lon_out, lat_out, verbose, &
                               interp_method, num_nbrs, max_dist, src_modulo,    &
                               mask_in, mask_out, is_latlon_in = src_is_latlon)

       call horiz_interp ( Interp, data_in, data_out, verbose,   &
                           missing_value=missing_value, missing_permit=missing_permit )
    end if

    call horiz_interp_del ( Interp )

!-----------------------------------------------------------------------

 end subroutine horiz_interp_solo_1d_dst_r8

!#######################################################################

!> Overloaded version of interface horiz_interp_solo_2
 subroutine horiz_interp_solo_old_r8 (data_in, wb, sb, dx, dy,  &
                                   lon_out, lat_out, data_out,  &
                                   verbose, mask_in, mask_out)

!-----------------------------------------------------------------------
      real(r8_kind), intent(in),  dimension(:,:) :: data_in !< Global input data stored from west to east
                                        !! (1st dimension), south to north (2nd dimension)
      real(r8_kind), intent(in)                  :: wb !< Longitude (radians) that correspond to western-most
                                              !! boundary of grid box j=1 in array data_in
      real(r8_kind), intent(in)                  :: sb !< Latitude (radians) that correspond to western-most
                                              !! boundary of grid box j=1 in array data_in
      real(r8_kind), intent(in)                  :: dx !< Grid spacing (in radians) for the longitude axis
                                              !! (first dimension) for the input data
      real(r8_kind), intent(in)                  :: dy !< Grid spacing (in radians) for the latitude axis
                                              !! (first dimension) for the input data
      real(r8_kind), intent(in),  dimension(:)   :: lon_out !< The longitude edges (in radians) for output
                                        !! data grid boxes. The values are for adjacent grid boxes
                                        !! and must increase in value. If there are MLON grid boxes
                                        !! there must be MLON+1 edge values
      real(r8_kind), intent(in),  dimension(:)   :: lat_out !< The latitude edges (in radians) for output
                                        !! data grid boxes. The values are for adjacent grid boxes
                                        !! and may increase or decrease in value. If there are NLAT
                                        !! grid boxes there must be NLAT+1 edge values
      real(r8_kind), intent(out), dimension(:,:) :: data_out !< Output data on the output grid defined by grid box
   integer, intent(in),                   optional :: verbose
      real(r8_kind), intent(in),   dimension(:,:), optional :: mask_in
      real(r8_kind), intent(out),  dimension(:,:), optional :: mask_out
!-----------------------------------------------------------------------
     real(r8_kind), dimension(size(data_in,1)+1)  :: blon_in
     real(r8_kind), dimension(size(data_in,2)+1)  :: blat_in
     integer :: i, j, nlon_in, nlat_in
     real(r8_kind)    :: tpi
     integer, parameter    :: kindl = r8_kind !< real size at compile time
!-----------------------------------------------------------------------
   call horiz_interp_init

   tpi = 2.0_kindl * real(pi, r8_kind)
   nlon_in = size(data_in,1)
   nlat_in = size(data_in,2)

   do i = 1, nlon_in+1
      blon_in(i) = wb + real(i-1, r8_kind)*dx
   enddo
      if (abs(blon_in(nlon_in+1)-blon_in(1)-tpi) < epsilon(blon_in)) &
              blon_in(nlon_in+1)=blon_in(1)+tpi

   do j = 2, nlat_in
      blat_in(j) = sb + real(j-1, r8_kind)*dy
   enddo
      blat_in(1)         = -0.5_kindl * real(pi, r8_kind)
      blat_in(nlat_in+1) =  0.5_kindl * real(pi, r8_kind)


   call horiz_interp_solo_1d (data_in, blon_in, blat_in,    &
                              lon_out, lat_out, data_out,   &
                              verbose, mask_in, mask_out    )

!-----------------------------------------------------------------------

 end subroutine horiz_interp_solo_old_r8

!#######################################################################


 !####################################################################
 function is_lat_lon_r8(lon, lat)
    real(r8_kind), dimension(:,:), intent(in) :: lon, lat
    logical                          :: is_lat_lon_r8
    integer                          :: i, j, nlon, nlat, num

    is_lat_lon_r8 = .true.
    nlon = size(lon,1)
    nlat = size(lon,2)
    LOOP_LAT: do j = 1, nlat
       do i = 2, nlon
          if(lat(i,j) .NE. lat(1,j)) then
             is_lat_lon_r8 = .false.
             exit LOOP_LAT
          end if
       end do
    end do LOOP_LAT

    if(is_lat_lon_r8) then
       LOOP_LON: do i = 1, nlon
          do j = 2, nlat
             if(lon(i,j) .NE. lon(i,1)) then
                is_lat_lon_r8 = .false.
                exit LOOP_LON
             end if
          end do
       end do LOOP_LON
    end if

    num = 0
    if(is_lat_lon_r8) num = 1
    call mpp_min(num)
    if(num == 1) then
       is_lat_lon_r8 = .true.
    else
       is_lat_lon_r8 = .false.
    end if

    return
 end function is_lat_lon_r8

 !> Subroutine for reading a weight file and use it to fill in the horiz interp type
!! for the bilinear interpolation method.
 subroutine horiz_interp_read_weights_r8(Interp, weight_filename, lon_out, lat_out, lon_in, lat_in, &
                                       weight_file_source, interp_method, isw, iew, jsw, jew, nglon, nglat)
   type(horiz_interp_type), intent(inout) :: Interp             !< Horiz interp time to fill
   character(len=*),        intent(in)    :: weight_filename    !< Name of the weight file
   real(r8_kind),      intent(in)    :: lat_out(:,:)       !< Output (model) latitude
   real(r8_kind),      intent(in)    :: lon_out(:,:)       !< Output (model) longitude
   real(r8_kind),      intent(in)    :: lat_in(:)          !< Input (data) latitude
   real(r8_kind),      intent(in)    :: lon_in(:)          !< Input (data) longitude
   character(len=*),        intent(in)    :: weight_file_source !< Source of the weight file
   character(len=*),        intent(in)    :: interp_method      !< The interp method to use
   integer,                 intent(in)    :: isw, iew, jsw, jew !< Starting and ending indices of the compute domain
   integer,                 intent(in)    :: nglon              !< Number of longitudes in the global domain
   integer,                 intent(in)    :: nglat              !< Number of latitudes in the globl domain

   integer                           :: i, j                 !< For do loops
   integer                           :: nlon_in              !< Number of longitude in the data
   integer                           :: nlat_in              !< Number of latitude in the data grid
   real(r8_kind), allocatable   :: lon_src_1d(:)        !< Center points of the longitude data grid
   real(r8_kind), allocatable   :: lat_src_1d(:)        !< Center points of the lattiude data grid
   integer, parameter                :: kindl = r8_kind !< real kind size currently compiling

   select case (trim(interp_method))
   case ("bilinear")
     !! This is to reproduce the behavior in horiz_interp_new
     !! The subroutine assumes that the data grid (lon_in, lat_in) are
     !! the edges and not the centers.
     !! Data_override passes in the edges, which are calculated using the axis_edges subroutine
     nlon_in  = size(lon_in(:))-1;  nlat_in  = size(lat_in(:))-1
     allocate(lon_src_1d(nlon_in), lat_src_1d(nlat_in))
     do i = 1, nlon_in
       lon_src_1d(i) = (lon_in(i) + lon_in(i+1)) * 0.5_kindl
     enddo
     do j = 1, nlat_in
       lat_src_1d(j) = (lat_in(j) + lat_in(j+1)) * 0.5_kindl
     enddo

     call horiz_interp_read_weights_bilinear(Interp, weight_filename, lon_out, lat_out, &
      lon_src_1d, lat_src_1d, weight_file_source, interp_method, &
      isw, iew, jsw, jew, nglon, nglat)
     deallocate(lon_src_1d,lat_src_1d)
   case default
     call mpp_error(FATAL, "Reading weight from file is not supported for the "//&
      trim(interp_method)//" method. It is currently only supported for bilinear")
   end select
 end subroutine horiz_interp_read_weights_r8
# 64 "horiz_interp/include/horiz_interp_r8.fh" 2
# 319 "horiz_interp/horiz_interp.F90" 2

end module horiz_interp_mod
!> @}
! close documentation grouping
