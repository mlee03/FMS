# 1 "topography/gaussian_topog.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "topography/gaussian_topog.F90"
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
!> @defgroup gaussian_topog_mod gaussian_topog_mod
!! @ingroup topography
!! @{
!! @brief Routines for creating Gaussian-shaped land surface topography
!! for latitude-longitude grids.
!! @author Bruce Wyman
!!
!! Interfaces generate simple Gaussian-shaped mountains from
!! parameters specified by either argument list or namelist input.
!! The mountain shapes are controlled by the height, half-width,
!! and ridge-width parameters.

module gaussian_topog_mod

use  fms_mod, only: check_nml_error,                 &
                    stdlog, write_version_number,    &
                    mpp_pe, mpp_root_pe,             &
                    error_mesg, FATAL

use constants_mod, only: pi

use mpp_mod,       only: input_nml_file
use platform_mod,  only: r4_kind, r8_kind

implicit none
private

public :: gaussian_topog_init, get_gaussian_topog

interface gaussian_topog_init
   module procedure gaussian_topog_init_r4
   module procedure gaussian_topog_init_r8
end interface gaussian_topog_init

interface get_gaussian_topog
    module procedure get_gaussian_topog_r4, get_gaussian_topog_r8
end interface get_gaussian_topog

!! Namelist information for gaussian_topog_nml
!!
!!     The variables in this namelist are only used when routine
!!     <TT>gaussian_topog_init</TT> is called.  The namelist variables
!!     are dimensioned (by 10), so that multiple mountains can be generated.
!!
!!     Internal parameter mxmtns = 10. By default no mountains are generated.
!!    </DATA>
!!
!!     NAMELIST FOR GENERATING GAUSSIAN MOUNTAINS
!!
!!  * multiple mountains can be generated
!!  * the final mountains are the sum of all
!!
!!       height = height in meters
!!       olon, olat = longitude,latitude origin              (degrees)
!!       rlon, rlat = longitude,latitude half-width of ridge (degrees)
!!       wlon, wlat = longitude,latitude half-width of tail  (degrees)
!!
!!       Note: For the standard gaussian mountain
!!             set rlon = rlat = 0 .
!!
!! <PRE>
!!
!!       height -->   ___________________________
!!                   /                           !!                  /              |              !!    gaussian     /               |               !!      sides --> /                |                !!               /               olon                !!         _____/                olat                 \______





!!
!!              |    |             |
!!              |<-->|<----------->|
!!              |wlon|    rlon     |
!!               wlat     rlat
!!
   integer, parameter :: maxmts = 10

   real(kind=r8_kind), dimension(maxmts) :: height = 0.0_r8_kind !< height in meters of the gaussian mountiains
   real(kind=r8_kind), dimension(maxmts) ::  olon  = 0.0_r8_kind !< longitude of mountain origins (degrees)
   real(kind=r8_kind), dimension(maxmts) ::  olat  = 0.0_r8_kind !< Latitude  of mountain origins (degrees)
   real(kind=r8_kind), dimension(maxmts) ::  wlon  = 0.0_r8_kind !< Longitude of half-width mountain trails (degrees)
   real(kind=r8_kind), dimension(maxmts) ::  wlat  = 0.0_r8_kind !< Latitude of half-width mountain trails (degrees)
   real(kind=r8_kind), dimension(maxmts) ::  rlon  = 0.0_r8_kind !< Longitude of half-width mountain ridges (degrees)
                                                                !! for "standard" gaussian mountain set, rlon/rlat = 0
   real(kind=r8_kind), dimension(maxmts) ::  rlat  = 0.0_r8_kind !< Latitude of half-width mountain ridges (degrees)
                                                                !! for "standard" gaussian mountain set, rlon/rlat = 0

   namelist /gaussian_topog_nml/ height, olon, olat, wlon, wlat, rlon, rlat
! </NAMELIST>

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
# 112 "topography/gaussian_topog.F90" 2

logical :: do_nml = .true.
logical :: module_is_initialized = .FALSE.

!-----------------------------------------------------------------------

contains

!#######################################################################

subroutine read_namelist

   integer :: iunit, ierr, io

!>  read namelist

   read (input_nml_file, gaussian_topog_nml, iostat=io)
   ierr = check_nml_error(io,'gaussian_topog_nml')

!>  write version and namelist to log file

   if (mpp_pe() == mpp_root_pe()) then
      iunit = stdlog()
      write (iunit, nml=gaussian_topog_nml)
   endif

   do_nml = .false.

end subroutine read_namelist

!#######################################################################


# 1 "topography/include/gaussian_topog_r4.fh" 1
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











# 1 "topography/include/gaussian_topog.inc" 1
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

!#######################################################################
!> @brief Returns a simple surface height field that consists of a single
!! Gaussian-shaped mountain.
!!
!> Returns a surface height field that consists
!! of the sum of one or more Gaussian-shaped mountains.
!!
subroutine gaussian_topog_init_r4 ( lon, lat, zsurf )

real(kind=r4_kind), intent(in)  :: lon(:) !< The mean grid box longitude in radians
real(kind=r4_kind), intent(in)  :: lat(:) !< The mean grid box latitude in radians
real(kind=r4_kind), intent(out) :: zsurf(:,:) !< The surface height (meters). Size must be size(lon) by size(lat)

integer :: n
integer, parameter :: lkind=r4_kind !local r4_kind kind

  if (.not.module_is_initialized) then
     call write_version_number("GAUSSIAN_TOPOG_MOD", version)
  endif

  if(any(shape(zsurf) /= (/size(lon(:)),size(lat(:))/))) then
    call error_mesg ('get_gaussian_topog in topography_mod', &
     'shape(zsurf) is not equal to (/size(lon),size(lat)/)', FATAL)
  endif

  if (do_nml) call read_namelist

! compute sum of all non-zero mountains
  zsurf(:,:) = 0.0_lkind
  do n = 1, maxmts
    if ( height(n) == 0.0_r8_kind ) cycle
    zsurf = zsurf + get_gaussian_topog ( lon, lat, real(height(n),lkind), &
                real(olon(n),lkind), real(olat(n),lkind), real(wlon(n),lkind), &
                real(wlat(n),lkind), real(rlon(n),lkind), real(rlat(n),lkind))
  enddo
 module_is_initialized = .TRUE.

end subroutine gaussian_topog_init_r4
!#######################################################################
!> The height, position, width, and elongation of the mountain
!! is controlled by optional arguments.
!! @param real lon The mean grid box longitude in radians.
!! @param real lat The mean grid box latitude in radians.
!! @param real height Maximum surface height in meters.
!! @param real olond, olatd Position/origin of mountain in degrees longitude and latitude.
!! This is the location of the maximum height.
!! @param real wlond, wlatd Gaussian half-width of mountain in degrees longitude and latitude.
!! @param real rlond, rlatd Ridge half-width of mountain in degrees longitude and latitude.
!! This is the elongation of the maximum height.
!! @param real zsurf The surface height (in meters).
!! The size of the returned field is size(lon) by size(lat).
!!   </OUT>
!!
!! @throws FATAL shape(zsurf) is not equal to (/size(lon),size(lat)/)
!!     Check the input grid size and output field size.
!!     The input grid is defined at the midpoint of grid boxes.
!!
!! @note
!!     Mountains do not wrap around the poles.
!
!! <br>Example usage:
!! @code{.F90} zsurf = <B>get_gaussian_topog</B> ( lon, lat, height
!!                    [, olond, olatd, wlond, wlatd, rlond, rlatd ] )@endcode
!> Returns a land surface topography that consists of a "set" of
!! simple Gaussian-shaped mountains.  The height, position,
!! width, and elongation of the mountains can be controlled
!! by variables in the namelist.
function get_gaussian_topog_r4(lon, lat, height,                         &
                            olond, olatd, wlond, wlatd, rlond, rlatd ) &
                     result (zsurf )

real(kind=r4_kind), intent(in)  :: lon(:), lat(:)
real(kind=r4_kind), intent(in)  :: height
real(kind=r4_kind), intent(in), optional :: olond, olatd, wlond, wlatd, rlond, rlatd
real(kind=r4_kind) :: zsurf(size(lon,1),size(lat,1))

integer :: i, j
real(kind=r4_kind)    :: olon, olat, wlon, wlat, rlon, rlat
real(kind=r4_kind)    :: tpi, dtr, dx, dy, xx, yy
integer, parameter     :: lkind = r4_kind

  if (do_nml) call read_namelist

! no need to compute mountain if height=0
  if ( height == 0.0_lkind) then
       zsurf(:,:) = 0.0_lkind
       return
  endif

  tpi = 2.0_lkind*real(pi, r4_kind)
  dtr = tpi/360.0_lkind

! defaults and convert degrees to radians (dtr)
  olon = 90.0_r8_kind*real(dtr, r8_kind);  if (present(olond)) olon=real(olond*dtr,r8_kind)
  olat = 45.0_r8_kind*real(dtr, r8_kind);  if (present(olatd)) olat=real(olatd*dtr,r8_kind)
  wlon = 15.0_r8_kind*real(dtr, r8_kind);  if (present(wlond)) wlon=real(wlond*dtr,r8_kind)
  wlat = 15.0_r8_kind*real(dtr, r8_kind);  if (present(wlatd)) wlat=real(wlatd*dtr,r8_kind)
  rlon = 0.0_r8_kind     ;  if (present(rlond)) rlon=real(rlond*dtr,r8_kind)
  rlat = 0.0_r8_kind     ;  if (present(rlatd)) rlat=real(rlatd*dtr,r8_kind)

! compute gaussian-shaped mountain
    do j=1,size(lat(:))
      dy = abs(lat(j) - real(olat,lkind))   ! dist from y origin
      yy = max(0.0_lkind, dy-real(rlat,lkind))/real(wlat,lkind)
      do i=1,size(lon(:))
        dx = abs(lon(i) - real(olon,lkind)) ! dist from x origin
        dx = min(dx, abs(dx-tpi))  ! To ensure that: -pi <= dx <= pi
        xx = max(0.0_lkind, dx-real(rlon,lkind))/real(wlon,lkind)
        zsurf(i,j) = real(height,lkind)*exp(-xx**2 - yy**2)
      enddo
    enddo

end function get_gaussian_topog_r4
# 29 "topography/include/gaussian_topog_r4.fh" 2
# 145 "topography/gaussian_topog.F90" 2

# 1 "topography/include/gaussian_topog_r8.fh" 1
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










# 1 "topography/include/gaussian_topog.inc" 1
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

!#######################################################################
!> @brief Returns a simple surface height field that consists of a single
!! Gaussian-shaped mountain.
!!
!> Returns a surface height field that consists
!! of the sum of one or more Gaussian-shaped mountains.
!!
subroutine gaussian_topog_init_r8 ( lon, lat, zsurf )

real(kind=r8_kind), intent(in)  :: lon(:) !< The mean grid box longitude in radians
real(kind=r8_kind), intent(in)  :: lat(:) !< The mean grid box latitude in radians
real(kind=r8_kind), intent(out) :: zsurf(:,:) !< The surface height (meters). Size must be size(lon) by size(lat)

integer :: n
integer, parameter :: lkind=r8_kind !local r8_kind kind

  if (.not.module_is_initialized) then
     call write_version_number("GAUSSIAN_TOPOG_MOD", version)
  endif

  if(any(shape(zsurf) /= (/size(lon(:)),size(lat(:))/))) then
    call error_mesg ('get_gaussian_topog in topography_mod', &
     'shape(zsurf) is not equal to (/size(lon),size(lat)/)', FATAL)
  endif

  if (do_nml) call read_namelist

! compute sum of all non-zero mountains
  zsurf(:,:) = 0.0_lkind
  do n = 1, maxmts
    if ( height(n) == 0.0_r8_kind ) cycle
    zsurf = zsurf + get_gaussian_topog ( lon, lat, real(height(n),lkind), &
                real(olon(n),lkind), real(olat(n),lkind), real(wlon(n),lkind), &
                real(wlat(n),lkind), real(rlon(n),lkind), real(rlat(n),lkind))
  enddo
 module_is_initialized = .TRUE.

end subroutine gaussian_topog_init_r8
!#######################################################################
!> The height, position, width, and elongation of the mountain
!! is controlled by optional arguments.
!! @param real lon The mean grid box longitude in radians.
!! @param real lat The mean grid box latitude in radians.
!! @param real height Maximum surface height in meters.
!! @param real olond, olatd Position/origin of mountain in degrees longitude and latitude.
!! This is the location of the maximum height.
!! @param real wlond, wlatd Gaussian half-width of mountain in degrees longitude and latitude.
!! @param real rlond, rlatd Ridge half-width of mountain in degrees longitude and latitude.
!! This is the elongation of the maximum height.
!! @param real zsurf The surface height (in meters).
!! The size of the returned field is size(lon) by size(lat).
!!   </OUT>
!!
!! @throws FATAL shape(zsurf) is not equal to (/size(lon),size(lat)/)
!!     Check the input grid size and output field size.
!!     The input grid is defined at the midpoint of grid boxes.
!!
!! @note
!!     Mountains do not wrap around the poles.
!
!! <br>Example usage:
!! @code{.F90} zsurf = <B>get_gaussian_topog</B> ( lon, lat, height
!!                    [, olond, olatd, wlond, wlatd, rlond, rlatd ] )@endcode
!> Returns a land surface topography that consists of a "set" of
!! simple Gaussian-shaped mountains.  The height, position,
!! width, and elongation of the mountains can be controlled
!! by variables in the namelist.
function get_gaussian_topog_r8(lon, lat, height,                         &
                            olond, olatd, wlond, wlatd, rlond, rlatd ) &
                     result (zsurf )

real(kind=r8_kind), intent(in)  :: lon(:), lat(:)
real(kind=r8_kind), intent(in)  :: height
real(kind=r8_kind), intent(in), optional :: olond, olatd, wlond, wlatd, rlond, rlatd
real(kind=r8_kind) :: zsurf(size(lon,1),size(lat,1))

integer :: i, j
real(kind=r8_kind)    :: olon, olat, wlon, wlat, rlon, rlat
real(kind=r8_kind)    :: tpi, dtr, dx, dy, xx, yy
integer, parameter     :: lkind = r8_kind

  if (do_nml) call read_namelist

! no need to compute mountain if height=0
  if ( height == 0.0_lkind) then
       zsurf(:,:) = 0.0_lkind
       return
  endif

  tpi = 2.0_lkind*real(pi, r8_kind)
  dtr = tpi/360.0_lkind

! defaults and convert degrees to radians (dtr)
  olon = 90.0_r8_kind*real(dtr, r8_kind);  if (present(olond)) olon=real(olond*dtr,r8_kind)
  olat = 45.0_r8_kind*real(dtr, r8_kind);  if (present(olatd)) olat=real(olatd*dtr,r8_kind)
  wlon = 15.0_r8_kind*real(dtr, r8_kind);  if (present(wlond)) wlon=real(wlond*dtr,r8_kind)
  wlat = 15.0_r8_kind*real(dtr, r8_kind);  if (present(wlatd)) wlat=real(wlatd*dtr,r8_kind)
  rlon = 0.0_r8_kind     ;  if (present(rlond)) rlon=real(rlond*dtr,r8_kind)
  rlat = 0.0_r8_kind     ;  if (present(rlatd)) rlat=real(rlatd*dtr,r8_kind)

! compute gaussian-shaped mountain
    do j=1,size(lat(:))
      dy = abs(lat(j) - real(olat,lkind))   ! dist from y origin
      yy = max(0.0_lkind, dy-real(rlat,lkind))/real(wlat,lkind)
      do i=1,size(lon(:))
        dx = abs(lon(i) - real(olon,lkind)) ! dist from x origin
        dx = min(dx, abs(dx-tpi))  ! To ensure that: -pi <= dx <= pi
        xx = max(0.0_lkind, dx-real(rlon,lkind))/real(wlon,lkind)
        zsurf(i,j) = real(height,lkind)*exp(-xx**2 - yy**2)
      enddo
    enddo

end function get_gaussian_topog_r8
# 28 "topography/include/gaussian_topog_r8.fh" 2
# 146 "topography/gaussian_topog.F90" 2

end module gaussian_topog_mod

!> @}
! close documentation grouping
