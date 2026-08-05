# 1 "horiz_interp/horiz_interp_spherical.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "horiz_interp/horiz_interp_spherical.F90"
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
!> @defgroup horiz_interp_spherical_mod horiz_interp_spherical_mod
!! @ingroup horiz_interp
!! @{
!! @brief Performs spatial interpolation between grids using inverse-distance-weighted scheme.
!! This module can interpolate data from rectangular/tripolar grid
!! to rectangular/tripolar grid. The interpolation scheme is inverse-distance-weighted
!! scheme.    There is an optional mask field for missing input data.
!! An optional output mask field may be used in conjunction with
!! the input mask to show where output data exists.

module horiz_interp_spherical_mod

  use platform_mod,          only : r4_kind, r8_kind
  use mpp_mod,               only : mpp_error, FATAL, WARNING, stdout
  use mpp_mod,               only : mpp_root_pe, mpp_pe
  use mpp_mod,               only : input_nml_file
  use fms_mod,               only : write_version_number
  use fms_mod,               only : check_nml_error
  use constants_mod,         only : pi
  use horiz_interp_type_mod, only : horiz_interp_type, stats, SPHERICAL

  implicit none
  private

  interface horiz_interp_spherical
    module procedure horiz_interp_spherical_r4
    module procedure horiz_interp_spherical_r8
  end interface

  interface horiz_interp_spherical_new
    module procedure horiz_interp_spherical_new_r4
    module procedure horiz_interp_spherical_new_r8
  end interface

  interface horiz_interp_spherical_wght
    module procedure horiz_interp_spherical_wght_r4
    module procedure horiz_interp_spherical_wght_r8
  end interface

  public :: horiz_interp_spherical_new, horiz_interp_spherical, horiz_interp_spherical_del
  public :: horiz_interp_spherical_init, horiz_interp_spherical_wght

 !> private helper routines
  interface full_search
    module procedure full_search_r4
    module procedure full_search_r8
  end interface

  interface radial_search
    module procedure radial_search_r4
    module procedure radial_search_r8
  end interface

  interface spherical_distance
    module procedure spherical_distance_r4
    module procedure spherical_distance_r8
  end interface

  integer, parameter :: max_neighbors = 400
  real(R8_KIND),    parameter :: max_dist_default = 0.1_r8_kind  ! radians
  integer, parameter :: num_nbrs_default = 4
  real(R8_KIND),    parameter :: large=1.e20_r8_kind
  real(R8_KIND),    parameter :: epsln=1.e-10_r8_kind

  integer            :: pe, root_pe


  character(len=32) :: search_method = "radial_search" !< Namelist variable to indicate the searching
                    !! method to find the
                    !! nearest neighbor points. Its value can be "radial_search" and "full_search",
                    !! with default value "radial_search". when search_method is "radial_search",
                    !! the search may be not quite accurate for some cases. Normally the search will
                    !! be ok if you chose suitable max_dist. When search_method is "full_search",
                    !! it will be always accurate, but will be slower comparing to "radial_search".
                    !! Normally these two search algorithm will produce same results other than
                    !! order of operation. "radial_search" are recommended to use. The purpose to
                    !! add "full_search" is in case you think you interpolation results is
                    !! not right, you have other option to verify.

!or "full_search"
  namelist /horiz_interp_spherical_nml/ search_method

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
# 103 "horiz_interp/horiz_interp_spherical.F90" 2
  logical            :: module_is_initialized = .FALSE.

contains

  !#######################################################################

  !> Initializes module and writes version number to logfile.out
  subroutine horiz_interp_spherical_init
    integer :: ierr, io


    if(module_is_initialized) return
    call write_version_number("horiz_interp_spherical_mod", version)
    read (input_nml_file, horiz_interp_spherical_nml, iostat=io)
    ierr = check_nml_error(io,'horiz_interp_spherical_nml')

     module_is_initialized = .true.

  end subroutine horiz_interp_spherical_init

  !#######################################################################

  !> Deallocates memory used by "HI_KIND_TYPE" variables.
  !! Must be called before reinitializing with horiz_interp_spherical_new.
  subroutine horiz_interp_spherical_del( Interp )

    type (horiz_interp_type), intent(inout) :: Interp !< A derived-type variable returned by previous
                                           !! call to horiz_interp_spherical_new. The input variable
                                           !! must have allocated arrays. The returned variable will
                                           !! contain deallocated arrays.

    if(Interp%horizInterpReals4_type%is_allocated) then
      if(allocated(Interp%horizInterpReals4_type%src_dist)) deallocate(Interp%horizInterpReals4_type%src_dist)
    else if (Interp%horizInterpReals8_type%is_allocated) then
      if(allocated(Interp%horizInterpReals8_type%src_dist)) deallocate(Interp%horizInterpReals8_type%src_dist)
    endif
    if(allocated(Interp%num_found)) deallocate(Interp%num_found)
    if(allocated(Interp%i_lon))     deallocate(Interp%i_lon)
    if(allocated(Interp%j_lat))     deallocate(Interp%j_lat)

    Interp%horizInterpReals4_type%is_allocated = .false.
    Interp%horizInterpReals8_type%is_allocated = .false.

  end subroutine horiz_interp_spherical_del

  !#######################################################################


# 1 "horiz_interp/include/horiz_interp_spherical_r4.fh" 1
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




























# 1 "horiz_interp/include/horiz_interp_spherical.inc" 1
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
  !> Initialization routine.
  !!
  !> Allocates space and initializes a derived-type variable
  !! that contains pre-computed interpolation indices and weights.
  subroutine horiz_interp_spherical_new_r4(Interp, lon_in,lat_in,lon_out,lat_out, &
       num_nbrs, max_dist, src_modulo)

    type(horiz_interp_type), intent(inout) :: Interp !< A derived type variable containing indices
                                          !! and weights for subsequent interpolations. To
                                          !! reinitialize for different grid-to-grid interpolation
                                          !! @ref horiz_interp_del must be used first.
    real(r4_kind), intent(in),  dimension(:,:)      :: lon_in !< Latitude (radians) for source data grid
    real(r4_kind), intent(in),  dimension(:,:)      :: lat_in !< Longitude (radians) for source data grid
    real(r4_kind), intent(in),  dimension(:,:)      :: lon_out !< Longitude (radians) for output data grid
    real(r4_kind), intent(in),  dimension(:,:)      :: lat_out !< Latitude (radians) for output data grid
    logical, intent(in),          optional :: src_modulo !< indicates if the boundary condition
                                          !! along zonal boundary is cyclic or not. Cyclic when true
    integer, intent(in),        optional   :: num_nbrs !< Number of nearest neighbors for regridding
                                           !! When number of neighbors within the radius max_dist
                                           !! is less than num_nbrs, All the neighbors will be used
                                           !! to interpolate onto destination grid. when number of
                                           !! neighbors within the radius max_dist is greater than
                                           !! num_nbrs, at least "num_nbrs"
  !     neighbors will be used to remap onto destination grid
    real(r4_kind), optional,             intent(in) :: max_dist !< Maximum region of influence around
                                                       !! destination grid points

    !------local variables ---------------------------------------
    integer :: i, j, n
    integer :: map_dst_xsize, map_dst_ysize, map_src_xsize, map_src_ysize
    integer :: map_src_size, num_neighbors
    real(r4_kind)    :: max_src_dist, tpi, hpi
    logical :: src_is_modulo
    real(r4_kind)    :: min_theta_dst, max_theta_dst, min_phi_dst, max_phi_dst
    real(r4_kind)    :: min_theta_src, max_theta_src, min_phi_src, max_phi_src
    integer               :: map_src_add(size(lon_out,1),size(lon_out,2),max_neighbors)
    real(r4_kind)    :: map_src_dist(size(lon_out,1),size(lon_out,2),max_neighbors)
    integer               :: num_found(size(lon_out,1),size(lon_out,2))
    integer               :: ilon(max_neighbors), jlat(max_neighbors)
    real(r4_kind), dimension(size(lon_out,1),size(lon_out,2)) :: theta_dst, phi_dst
    real(r4_kind), dimension(size(lon_in,1)*size(lon_in,2))   :: theta_src, phi_src
    integer, parameter                                             :: kindl = r4_kind

    !--------------------------------------------------------------

    pe      = mpp_pe()
    root_pe = mpp_root_pe()

    tpi = 2.0_kindl*real(PI, r4_kind); hpi = 0.5_kindl*real(PI,r4_kind)

    num_neighbors = num_nbrs_default
    if(present(num_nbrs)) num_neighbors = num_nbrs
    if (num_neighbors <= 0) call mpp_error(FATAL,'horiz_interp_spherical_mod: num_neighbors must be > 0')

    max_src_dist = real(max_dist_default, r4_kind)
    if (PRESENT(max_dist)) max_src_dist = max_dist
    Interp%horizInterpReals4_type%max_src_dist = max_src_dist

    src_is_modulo = .true.
    if (PRESENT(src_modulo)) src_is_modulo = src_modulo

    !--- check the grid size comformable
    map_dst_xsize=size(lon_out,1);map_dst_ysize=size(lon_out,2)
    map_src_xsize=size(lon_in,1); map_src_ysize=size(lon_in,2)
    map_src_size = map_src_xsize*map_src_ysize

    if (map_dst_xsize /= size(lat_out,1) .or. map_dst_ysize /= size(lat_out,2)) &
         call mpp_error(FATAL,'horiz_interp_spherical_mod: destination grids not conformable')
    if (map_src_xsize /= size(lat_in,1) .or. map_src_ysize /= size(lat_in,2)) &
         call mpp_error(FATAL,'horiz_interp_spherical_mod: source grids not conformable')

    theta_src      = reshape(lon_in,(/map_src_size/))
    phi_src        = reshape(lat_in,(/map_src_size/))
    theta_dst(:,:) = lon_out(:,:)
    phi_dst(:,:)   = lat_out(:,:)

    min_theta_dst=real(tpi, r4_kind);max_theta_dst=0.0_kindl
    min_phi_dst=real(pi, r4_kind);max_phi_dst=real(-pi, r4_kind)
    min_theta_src=real(tpi, r4_kind);max_theta_src=0.0_kindl
    min_phi_src=real(pi, r4_kind);max_phi_src=real(-pi, r4_kind)

    where(theta_dst<0.0_kindl)  theta_dst = theta_dst+real(tpi,r4_kind)

    where(theta_dst>real(tpi,r4_kind))  theta_dst = theta_dst-real(tpi,r4_kind)
    where(theta_src<0.0_kindl)  theta_src =  theta_src+real(tpi,r4_kind)
    where(theta_src>real(tpi,r4_kind))  theta_src = theta_src-real(tpi,r4_kind)

    where(phi_dst < -hpi) phi_dst = -hpi
    where(phi_dst > hpi)  phi_dst =  hpi
    where(phi_src < -hpi) phi_src = -hpi
    where(phi_src > hpi)  phi_src =  hpi

    do j=1,map_dst_ysize
       do i=1,map_dst_xsize
          min_theta_dst = min(min_theta_dst,theta_dst(i,j))
          max_theta_dst = max(max_theta_dst,theta_dst(i,j))
          min_phi_dst = min(min_phi_dst,phi_dst(i,j))
          max_phi_dst = max(max_phi_dst,phi_dst(i,j))
       enddo
    enddo

    do i=1,map_src_size
       min_theta_src = min(min_theta_src,theta_src(i))
       max_theta_src = max(max_theta_src,theta_src(i))
       min_phi_src = min(min_phi_src,phi_src(i))
       max_phi_src = max(max_phi_src,phi_src(i))
    enddo

    if (min_phi_dst < min_phi_src) print *, '=> WARNING:  latitute of dest grid exceeds src'
    if (max_phi_dst > max_phi_src) print *, '=> WARNING:  latitute of dest grid exceeds src'
    ! when src is cyclic, no need to print out the following warning.
    if(.not. src_is_modulo) then
       if (min_theta_dst < min_theta_src) print *, '=> WARNING : longitude of dest grid exceeds src'
       if (max_theta_dst > max_theta_src) print *, '=> WARNING : longitude of dest grid exceeds src'
    endif

    ! allocate memory to data type
    if(ALLOCATED(Interp%i_lon)) then
       if(size(Interp%i_lon,1) .NE. map_dst_xsize .OR.     &
          size(Interp%i_lon,2) .NE. map_dst_ysize )  call mpp_error(FATAL, &
          'horiz_interp_spherical_mod: size(Interp%i_lon(:),1) .NE. map_dst_xsize .OR. '// &
          'size(Interp%i_lon(:),2) .NE. map_dst_ysize')
    else
       allocate(Interp%i_lon(map_dst_xsize,map_dst_ysize,max_neighbors), &
            Interp%j_lat(map_dst_xsize,map_dst_ysize,max_neighbors),     &
            Interp%horizInterpReals4_type%src_dist(map_dst_xsize,map_dst_ysize,max_neighbors),  &
            Interp%num_found(map_dst_xsize,map_dst_ysize)       )
    endif

    map_src_add         = 0
    map_src_dist        = real(large, r4_kind)
    num_found           = 0

    !using radial_search to find the nearest points and corresponding distance.

    select case(trim(search_method))
    case ("radial_search") ! will be efficient, but may be not so accurate for some cases
       call radial_search(theta_src, phi_src, theta_dst, phi_dst, map_src_xsize, map_src_ysize, &
            map_src_add, map_src_dist, num_found, num_neighbors,max_src_dist,src_is_modulo)
    case ("full_search")   ! always accurate, but less efficient.
       call full_search(theta_src, phi_src, theta_dst, phi_dst, map_src_add, map_src_dist, &
            num_found, num_neighbors,max_src_dist )
    case default
       call mpp_error(FATAL,"HORIZ_INTERP_SPHERICAL_NEW_: nml search_method = "// &
                  trim(search_method)//" is not a valid namelist option")
    end select

    do j=1,map_dst_ysize
       do i=1,map_dst_xsize
          do n=1,num_found(i,j)
             if(map_src_add(i,j,n) == 0) then
                jlat(n) = 0; ilon(n) = 0
             else
                jlat(n) = map_src_add(i,j,n)/map_src_xsize + 1
                ilon(n) = map_src_add(i,j,n) - (jlat(n)-1)*map_src_xsize
                if(ilon(n) == 0) then
                   jlat(n) = jlat(n) - 1
                   ilon(n) = map_src_xsize
                endif
             endif
          enddo
          Interp%i_lon(i,j,:)         = ilon(:)
          Interp%j_lat(i,j,:)         = jlat(:)
          Interp%num_found(i,j)       = num_found(i,j)
          Interp%horizInterpReals4_type%src_dist(i,j,:)      = map_src_dist(i,j,:)
       enddo
    enddo

    Interp%nlon_src = map_src_xsize; Interp%nlat_src = map_src_ysize
    Interp%nlon_dst = map_dst_xsize; Interp%nlat_dst = map_dst_ysize
    Interp% horizInterpReals4_type % is_allocated = .true.
    Interp% interp_method = SPHERICAL

    return

  end subroutine horiz_interp_spherical_new_r4

  !#######################################################################

  !> Subroutine for performing the horizontal interpolation between two grids.
  !! horiz_interp_spherical_new_r4 must be called before calling this routine.
  subroutine horiz_interp_spherical_r4( Interp, data_in, data_out, verbose, mask_in, mask_out, missing_value)
    type(horiz_interp_type), intent(in) :: Interp !< A derived type variable containing indices
                                          !! and weights for subsequent interpolations. Returned
                                          !! by a previous call to horiz_interp_spherical_new_r4
    real(r4_kind), intent(in),  dimension(:,:)           :: data_in !< Input data on source grid
    real(r4_kind), intent(out), dimension(:,:)           :: data_out !< Output data on destination grid
    integer, intent(in),               optional :: verbose !< 0 = no output; 1 = min,max,means; 2 = most output
    real(r4_kind), intent(in), dimension(:,:),  optional :: mask_in !< Input mask, must be the same size as
                                             !! the input data. The real(r4_kind) value of mask_in must be
                                             !! in the range (0.,1.). Set mask_in=0.0 for data points
                                             !! that should not be used or have missing data
    real(r4_kind), intent(out), dimension(:,:), optional :: mask_out !< Output mask that specifies whether data
                                                                         !! was computed.
    real(r4_kind), intent(in),                  optional :: missing_value !< Used to indicate missing data

    !--- some local variables ----------------------------------------
    real(r4_kind), dimension(Interp%nlon_dst, Interp%nlat_dst,size(Interp%horizInterpReals4_type%src_dist,3)) :: wt
    real(r4_kind), dimension(Interp%nlon_src, Interp%nlat_src) :: mask_src
    real(r4_kind), dimension(Interp%nlon_dst, Interp%nlat_dst) :: mask_dst
    integer :: nlon_in, nlat_in, nlon_out, nlat_out, num_found
    integer :: m, n, i, j, k, miss_in, miss_out, i1, i2, j1, j2, iverbose
    real(r4_kind)    :: min_in, max_in, avg_in, min_out, max_out, avg_out, sum
    integer, parameter    :: kindl = r4_kind !< compiled real kind size
    !-----------------------------------------------------------------

    iverbose = 0;  if (present(verbose)) iverbose = verbose

    nlon_in  = Interp%nlon_src; nlat_in  = Interp%nlat_src
    nlon_out = Interp%nlon_dst; nlat_out = Interp%nlat_dst

    if(size(data_in,1) .ne. nlon_in .or. size(data_in,2) .ne. nlat_in ) &
         call mpp_error(FATAL,'horiz_interp_spherical_mod: size of input array incorrect')

    if(size(data_out,1) .ne. nlon_out .or. size(data_out,2) .ne. nlat_out ) &
         call mpp_error(FATAL,'horiz_interp_spherical_mod: size of output array incorrect')

    mask_src = 1.0_kindl; mask_dst = 1.0_kindl
    if(present(mask_in)) mask_src = mask_in

    do n=1,nlat_out
       do m=1,nlon_out
          ! neighbors are sorted nearest to farthest
          ! check nearest to see if it is a land point
          num_found = Interp%num_found(m,n)
          if(num_found == 0 ) then
             mask_dst(m,n) = 0.0_kindl
          else
             i1 = Interp%i_lon(m,n,1); j1 = Interp%j_lat(m,n,1)
             if (mask_src(i1,j1) .lt. 0.5_kindl ) then
                mask_dst(m,n) = 0.0_kindl
             endif

             if(num_found .gt. 1 ) then
                i2 = Interp%i_lon(m,n,2); j2 = Interp%j_lat(m,n,2)
                ! compare first 2 nearest neighbors -- if they are nearly
                ! equidistant then use this mask for robustness
                if(abs(Interp%horizInterpReals4_type%src_dist(m,n,2)-Interp%horizInterpReals4_type%src_dist(m,n,1)) .lt. &
                       real(epsln,r4_kind)) then
                   if((mask_src(i1,j1) .lt. 0.5_kindl ))  mask_dst(m,n) = 0.0_kindl
                endif
             endif

             sum=0.0_kindl
             do k=1, num_found
                if(mask_src(Interp%i_lon(m,n,k),Interp%j_lat(m,n,k)) .lt. 0.5_kindl ) then
                   wt(m,n,k) = 0.0_kindl
                else
                   if (Interp%horizInterpReals4_type%src_dist(m,n,k) <= real(epsln,r4_kind)) then
                      wt(m,n,k) = real(large, r4_kind)
                      sum = sum + real(large, r4_kind)
                   else if(Interp%horizInterpReals4_type%src_dist(m,n,k) <= Interp%horizInterpReals4_type%max_src_dist ) then
                      wt(m,n,k) = 1.0_kindl /Interp%horizInterpReals4_type%src_dist(m,n,k)
                      sum = sum+wt(m,n,k)
                   else
                      wt(m,n,k) = 0.0_kindl
                   endif
                endif
             enddo
             if (sum > real(epsln,r4_kind)) then
                do k = 1, num_found
                   wt(m,n,k) = wt(m,n,k)/sum
                enddo
             else
                mask_dst(m,n) = 0.0_kindl
             endif
          endif
       enddo
    enddo

    data_out = 0.0_kindl
    do n=1,nlat_out
       do m=1,nlon_out
          if(mask_dst(m,n) .gt. 0.5_kindl ) then
             do k=1, Interp%num_found(m,n)
                i = Interp%i_lon(m,n,k)
                j = Interp%j_lat(m,n,k)
                data_out(m,n) = data_out(m,n)+data_in(i,j)*wt(m,n,k)
             enddo
          else
             if(present(missing_value)) then
                data_out(m,n) = missing_value
             else
                data_out(m,n) = 0.0_kindl
             endif
          endif
       enddo
    enddo

    if(present(mask_out)) mask_out = mask_dst

    !***********************************************************************
    ! compute statistics: minimum, maximum, and mean
    !-----------------------------------------------------------------------

    if (iverbose > 0) then

       ! compute statistics of input data

       call stats (data_in, min_in, max_in, avg_in, miss_in, missing_value, mask=mask_src)

       ! compute statistics of output data
       call stats (data_out, min_out, max_out, avg_out, miss_out, missing_value, mask=mask_dst)

       !---- output statistics ----
       ! root_pe have the information of global mean, min and max
       if(pe == root_pe) then
          write (*,900)
          write (*,901)  min_in ,max_in, avg_in
          if (present(mask_in))  write (*,903)  miss_in
          write (*,902)  min_out,max_out,avg_out
          if (present(mask_out)) write (*,903)  miss_out
       endif
900    format (/,1x,10('-'),' output from horiz_interp ',10('-'))
901    format ('  input:  min=',f16.9,'  max=',f16.9,'  avg=',f22.15)
902    format (' output:  min=',f16.9,'  max=',f16.9,'  avg=',f22.15)
903    format ('          number of missing points = ',i6)

    endif

    return
  end subroutine horiz_interp_spherical_r4

  !#######################################################################
  !> This routine isn't used internally
  !! it's similar to the routine above, just gets the weights as an out variable
  subroutine horiz_interp_spherical_wght_r4( Interp, wt, verbose, mask_in, mask_out, missing_value)
    type (horiz_interp_type), intent(in)        :: Interp
    real(r4_kind), intent(out), dimension(:,:,:)         :: wt
    integer, intent(in),               optional :: verbose
    real(r4_kind), intent(in), dimension(:,:),  optional :: mask_in
    real(r4_kind), intent(inout), dimension(:,:), optional :: mask_out
    real(r4_kind), intent(in),                  optional :: missing_value

    !--- some local variables ----------------------------------------
    real(r4_kind), dimension(Interp%nlon_src, Interp%nlat_src) :: mask_src
    real(r4_kind), dimension(Interp%nlon_dst, Interp%nlat_dst) :: mask_dst
    integer :: nlon_in, nlat_in, nlon_out, nlat_out, num_found
    integer :: m, n, k, i1, i2, j1, j2, iverbose
    real(r4_kind)    :: sum
    integer, parameter    :: kindl = r4_kind
    !-----------------------------------------------------------------

    iverbose = 0;  if (present(verbose)) iverbose = verbose

    nlon_in  = Interp%nlon_src; nlat_in  = Interp%nlat_src
    nlon_out = Interp%nlon_dst; nlat_out = Interp%nlat_dst

    mask_src = 1.0_kindl ; mask_dst = 1.0_kindl
    if(present(mask_in)) mask_src = mask_in

    do n=1,nlat_out
       do m=1,nlon_out
          ! neighbors are sorted nearest to farthest
          ! check nearest to see if it is a land point
          num_found = Interp%num_found(m,n)

          if (num_found > num_nbrs_default) then
            if( iverbose .gt. 0) print *,'pe=',mpp_pe(),'num_found=',num_found
            num_found = num_nbrs_default
          end if

          if(num_found == 0 ) then
             mask_dst(m,n) = 0.0_kindl
          else
             i1 = Interp%i_lon(m,n,1); j1 = Interp%j_lat(m,n,1)
             if (mask_src(i1,j1) .lt. 0.5_kindl) then
                mask_dst(m,n) = 0.0_kindl
             endif

             if(num_found .gt. 1 ) then
                i2 = Interp%i_lon(m,n,2); j2 = Interp%j_lat(m,n,2)
                ! compare first 2 nearest neighbors -- if they are nearly
                ! equidistant then use this mask for robustness
                if(abs(Interp%horizInterpReals4_type%src_dist(m,n,2)-Interp%horizInterpReals4_type%src_dist(m,n,1)) .lt. &
                   real(epsln,r4_kind)) then
                   if((mask_src(i1,j1) .lt. 0.5_kindl ))  mask_dst(m,n) = 0.0_kindl
                endif
             endif

             sum=0.0_kindl
             do k=1, num_found
                if(mask_src(Interp%i_lon(m,n,k),Interp%j_lat(m,n,k)) .lt. 0.5_kindl ) then
                   wt(m,n,k) = 0.0_kindl
                else
                   if (Interp%horizInterpReals4_type%src_dist(m,n,k) <= real(epsln, r4_kind)) then
                      wt(m,n,k) = real(large, r4_kind)
                      sum = sum + real(large, r4_kind)
                   else if(Interp%horizInterpReals4_type%src_dist(m,n,k) <= Interp%horizInterpReals4_type%max_src_dist ) then
                      wt(m,n,k) = 1.0_kindl /Interp%horizInterpReals4_type%src_dist(m,n,k)
                      sum = sum+wt(m,n,k)
                   else
                      wt(m,n,k) = 0.0_kindl
                   endif
                endif
             enddo
             if (sum > real(epsln,r4_kind)) then
                do k = 1, num_found
                   wt(m,n,k) = wt(m,n,k)/sum
                enddo
             else
                mask_dst(m,n) = 0.0_kindl
             endif
          endif
       enddo
    enddo

    return
  end subroutine horiz_interp_spherical_wght_r4

  !#######################################################################

  subroutine radial_search_r4(theta_src,phi_src,theta_dst,phi_dst, map_src_xsize, map_src_ysize, &
       map_src_add, map_src_dist, num_found, num_neighbors,max_src_dist,src_is_modulo)
    real(r4_kind),    intent(in),    dimension(:)   :: theta_src, phi_src
    real(r4_kind),    intent(in),    dimension(:,:) :: theta_dst, phi_dst
    integer, intent(in)                    :: map_src_xsize, map_src_ysize
    integer, intent(out), dimension(:,:,:) :: map_src_add
    real(r4_kind),    intent(out), dimension(:,:,:) :: map_src_dist
    integer, intent(inout), dimension(:,:) :: num_found
    integer, intent(in)                    :: num_neighbors
    real(r4_kind),    intent(in)                    :: max_src_dist
    logical, intent(in)                    :: src_is_modulo

    !---------- local variables ----------------------------------------
    integer, parameter :: max_nbrs = 50
    integer :: i, j, jj, i0, j0, n, l,i_left, i_right
    integer :: map_dst_xsize, map_dst_ysize
    integer :: i_left1, i_left2, i_right1, i_right2
    integer :: map_src_size, step, step_size, bound, bound_start, bound_end
    logical :: continue_search, result, continue_radial_search
    real(r4_kind)    :: d, res
    !------------------------------------------------------------------
    map_dst_xsize=size(theta_dst,1);map_dst_ysize=size(theta_dst,2)
    map_src_size = map_src_xsize*map_src_ysize

    do j=1,map_dst_ysize
       do i=1,map_dst_xsize
          continue_search=.true.
          step = 1
          step_size = int( sqrt(real(map_src_size, kind=r4_kind )))
          do while (continue_search .and. step_size > 0)
             do while (step <= map_src_size .and. continue_search)
                ! count land points as nearest neighbors
                d = spherical_distance(theta_dst(i,j),phi_dst(i,j),theta_src(step),phi_src(step))
                if (d <= max_src_dist) then
                   result = update_dest_neighbors_r4(map_src_add(i,j,:),map_src_dist(i,j,:), &
                        step,d, num_found(i,j), num_neighbors )
                   if (result) then
                      n = 0
                      i0 = mod(step,map_src_xsize)

                      if (i0 == 0) i0 = map_src_xsize
                      res = real(step, r4_kind)/real(map_src_xsize, r4_kind)
                      j0 = ceiling(res)
                      continue_radial_search = .true.
                      do while (continue_radial_search)
                         continue_radial_search = .false.
                         n = n+1 ! radial counter
                         if(n > max_nbrs) exit
                         ! ************** left boundary *******************************
                         i_left = i0-n
                         if (i_left <= 0) then
                            if (src_is_modulo) then
                               i_left = map_src_xsize + i_left
                            else
                               i_left = 1
                            endif
                         endif

                         do l = 0, 2*n
                            jj = j0 - n - 1 + l
                            if( jj < 0) then
                               bound = ( 1 - jj )*map_src_xsize - i_left
                            else if ( jj >= map_src_ysize ) then
                               bound = ( 2*map_src_ysize - jj ) * map_src_xsize - i_left
                            else
                               bound = jj * map_src_xsize + i_left
                            endif

                            d = spherical_distance(theta_dst(i,j),phi_dst(i,j),theta_src(bound), &
                                                                 phi_src(bound))
                            if(d<=max_src_dist) then
                               result = update_dest_neighbors_r4(map_src_add(i,j,:),map_src_dist(i,j,:), &
                                    bound,d, num_found(i,j), num_neighbors)
                               if (result) continue_radial_search = .true.
                            endif
                         enddo

                         ! ***************************right boundary *******************************
                         i_right = i0+n
                         if (i_right > map_src_xsize) then
                            if (src_is_modulo) then
                               i_right = i_right - map_src_xsize
                            else
                               i_right = map_src_xsize
                            endif
                         endif

                         do l = 0, 2*n
                            jj = j0 - n - 1 + l
                            if( jj < 0) then
                               bound = ( 1 - jj )*map_src_xsize - i_right
                            else if ( jj >= map_src_ysize ) then
                               bound = ( 2*map_src_ysize - jj) * map_src_xsize - i_right

                            else
                               bound = jj * map_src_xsize + i_right
                            endif

                            d = spherical_distance(theta_dst(i,j),phi_dst(i,j),theta_src(bound), &
                                                                 phi_src(bound))
                            if(d<=max_src_dist) then
                               result = update_dest_neighbors_r4(map_src_add(i,j,:),map_src_dist(i,j,:), &
                                    bound,d, num_found(i,j), num_neighbors)
                               if (result) continue_radial_search = .true.
                            endif
                         enddo

                         ! ************************* bottom boundary **********************************
                         i_left2 = 0
                         if( i_left > i_right) then
                            i_left1 = 1
                            i_right1 = i_right
                            i_left2 = i_left
                            i_right2 = map_src_xsize
                         else
                            i_left1 = i_left
                            i_right1 = i_right
                         endif

                         jj = j0 - n - 1
                         if( jj < 0 ) then
                            bound_start = ( 1 - jj)*map_src_xsize - i_right1
                            bound_end   = ( 1 - jj)*map_src_xsize - i_left1
                         else
                            bound_start = jj * map_src_xsize + i_left1
                            bound_end = jj * map_src_xsize + i_right1
                         endif

                         bound = bound_start
                         do while (bound <= bound_end)
                            d = spherical_distance(theta_dst(i,j),phi_dst(i,j),theta_src(bound), &
                                                                 phi_src(bound))
                            if(d<=max_src_dist) then
                               result = update_dest_neighbors_r4(map_src_add(i,j,:),map_src_dist(i,j,:), &
                                    bound,d, num_found(i,j), num_neighbors)
                               if (result) continue_radial_search = .true.
                            endif
                            bound = bound + 1

                         enddo

                         if(i_left2 > 0 ) then
                            if( jj < 0 ) then
                               bound_start = ( 1 - jj)*map_src_xsize - i_right2
                               bound_end   = ( 1 - jj)*map_src_xsize - i_left2
                            else
                               bound_start = jj * map_src_xsize + i_left2
                               bound_end = jj * map_src_xsize + i_right2
                            endif

                            bound = bound_start
                            do while (bound <= bound_end)
                               d = spherical_distance(theta_dst(i,j),phi_dst(i,j),theta_src(bound), &
                                                                    phi_src(bound))
                               if(d<=max_src_dist) then
                                  result = update_dest_neighbors_r4(map_src_add(i,j,:),map_src_dist(i,j,:), &
                                       bound,d, num_found(i,j), num_neighbors)
                                  if (result) continue_radial_search = .true.
                               endif
                               bound = bound + 1
                            enddo
                         endif

                         ! ************************** top boundary ************************************
                         jj = j0 + n - 1
                         if( jj >= map_src_ysize) then
                            bound_start = ( 2*map_src_ysize - jj ) * map_src_xsize - i_right1
                            bound_end   = ( 2*map_src_ysize - jj ) * map_src_xsize - i_left1
                         else
                            bound_start = jj * map_src_xsize + i_left1
                            bound_end = jj * map_src_xsize + i_right1
                         endif

                         bound = bound_start
                         do while (bound <= bound_end)
                            d = spherical_distance(theta_dst(i,j),phi_dst(i,j),theta_src(bound), &
                                                                 phi_src(bound))
                            if(d<=max_src_dist) then
                               result = update_dest_neighbors_r4(map_src_add(i,j,:),map_src_dist(i,j,:), &
                                    bound,d, num_found(i,j), num_neighbors)
                               if (result) continue_radial_search = .true.
                            endif
                            bound = bound + 1
                         enddo

                         if(i_left2 > 0) then
                            if( jj >= map_src_ysize) then
                               bound_start = ( 2*map_src_ysize - jj ) * map_src_xsize - i_right2
                               bound_end   = ( 2*map_src_ysize - jj ) * map_src_xsize - i_left2
                            else
                               bound_start = jj * map_src_xsize + i_left2
                               bound_end = jj * map_src_xsize + i_right2
                            endif

                            bound = bound_start
                            do while (bound <= bound_end)
                               d = spherical_distance(theta_dst(i,j),phi_dst(i,j),theta_src(bound), &
                                                                    phi_src(bound))
                               if(d<=max_src_dist) then
                                  result = update_dest_neighbors_r4(map_src_add(i,j,:),map_src_dist(i,j,:), &
                                       bound,d, num_found(i,j), num_neighbors)
                                  if (result) continue_radial_search = .true.
                               endif
                               bound = bound + 1
                            enddo
                         endif

                      enddo
                      continue_search = .false. ! stop looking
                   endif
                endif
                step=step+step_size
             enddo ! search loop
             step = 1
             step_size = step_size/2
          enddo
       enddo
    enddo

    return

  end subroutine radial_search_r4


  !#####################################################################

  function update_dest_neighbors_r4(map_src_add, map_src_dist, src_add,d, num_found, min_nbrs)

    integer, intent(inout), dimension(:) :: map_src_add
    real(r4_kind), intent(inout),    dimension(:) :: map_src_dist
    integer, intent(in)                  :: src_add
    real(r4_kind), intent(in)                     :: d
    integer, intent(inout)               :: num_found
    integer, intent(in)                  :: min_nbrs

    logical :: update_dest_neighbors_r4, already_exist = .false.

    integer :: n,m

    update_dest_neighbors_r4 = .false.

    n = 0
    NLOOP : do while ( n .le. num_found )
       n = n + 1
       DIST_CHK : if (d .le. map_src_dist(n)) then
          do m=n,num_found
             if (src_add == map_src_add(m)) then
                already_exist = .true.
                exit NLOOP
             endif
          enddo
          if(num_found < max_neighbors) then
             num_found = num_found + 1
          else
             call mpp_error(FATAL,'UPDATE_DEST_NEIGHBORS_: '// &
                  'number of neighbor points found is greated than maxium neighbor points' )
          endif
          do m=num_found,n+1,-1
             map_src_add(m) = map_src_add(m-1)
             map_src_dist(m) = map_src_dist(m-1)
          enddo
          map_src_add(n) = src_add
          map_src_dist(n) = d
          update_dest_neighbors_r4 = .true.
          if( num_found > min_nbrs ) then
             if( map_src_dist(num_found) > map_src_dist(num_found-1) ) then
                num_found = num_found - 1
             endif
             if( map_src_dist(min_nbrs+1) > map_src_dist(min_nbrs) ) then
                num_found = min_nbrs
             endif
          endif
          exit NLOOP ! n loop
       endif DIST_CHK
    end do NLOOP
    if(already_exist) return

    if( .not. update_dest_neighbors_r4 ) then
       if( num_found < min_nbrs ) then
          num_found               = num_found + 1
          update_dest_neighbors_r4   = .true.
          map_src_add(num_found)  = src_add
          map_src_dist(num_found) = d
       endif
    endif


    return

  end function update_dest_neighbors_r4

  !########################################################################
!  function spherical_distance_r4(theta1,phi1,theta2,phi2)

!    real(r4_kind), intent(in) :: theta1, phi1, theta2, phi2
!    real(r4_kind) :: spherical_distance_r4

!    real(r4_kind) :: r1(3), r2(3), cross(3), s, dot, ang

    ! this is a simple, enough way to calculate distance on the sphere
    ! first, construct cartesian vectors r1 and r2
    ! then calculate the cross-product which is proportional to the area
    ! between the 2 vectors.  The angular distance is arcsin of the
    ! distancealong the sphere
    !
    ! theta is longitude and phi is latitude
    !


!    r1(1) = cos(theta1)*cos(phi1);r1(2)=sin(theta1)*cos(phi1);r1(3)=sin(phi1)
!    r2(1) = cos(theta2)*cos(phi2);r2(2)=sin(theta2)*cos(phi2);r2(3)=sin(phi2)

!    cross(1) = r1(2)*r2(3)-r1(3)*r2(2)
!    cross(2) = r1(3)*r2(1)-r1(1)*r2(3)
!    cross(3) = r1(1)*r2(2)-r1(2)*r2(1)

!    s = sqrt(cross(1)**2.+cross(2)**2.+cross(3)**2.)

!    s = min(s,real(1.0, r4_kind)-epsln)

!    dot = r1(1)*r2(1) + r1(2)*r2(2) + r1(3)*r2(3)

!    if (dot > 0) then
!       ang = asin(s)
!    else if (dot < 0) then
!       ang = pi + asin(s)  !?  original is pi - asin(s)
!    else
!       ang = pi/real(2., r4_kind)
!    endif

!    spherical_distance_r4 = abs(ang) ! in radians

!    return

!  end function spherical_distance_r4
  ! The great cycle distance
  function spherical_distance_r4(theta1,phi1,theta2,phi2)

    real(r4_kind), intent(in) :: theta1, phi1, theta2, phi2
    real(r4_kind) :: spherical_distance_r4, dot
    integer, parameter :: kindl = r4_kind

    if(theta1 == theta2 .and. phi1 == phi2) then
        spherical_distance_r4 = 0.0_kindl
        return
    endif

    dot = cos(phi1)*cos(phi2)*cos(theta1-theta2) + sin(phi1)*sin(phi2)
    if(dot > 1.0_kindl)  dot = 1.0_kindl
    if(dot < real(-1.0_kindl, r4_kind)) dot = -1.0_kindl
    spherical_distance_r4 = acos(dot)

    return

  end function spherical_distance_r4


  !#######################################################################

  subroutine full_search_r4(theta_src,phi_src,theta_dst,phi_dst,map_src_add, map_src_dist,num_found, &
                         num_neighbors,max_src_dist)
    real(r4_kind),    intent(in),    dimension(:)   :: theta_src, phi_src
    real(r4_kind),    intent(in),    dimension(:,:) :: theta_dst, phi_dst
    integer, intent(out), dimension(:,:,:) :: map_src_add
    real(r4_kind),    intent(out), dimension(:,:,:) :: map_src_dist
    integer, intent(out), dimension(:,:)   :: num_found
    integer, intent(in)                    :: num_neighbors
    real(r4_kind), intent(in)                       :: max_src_dist

    integer :: i,j,map_src_size, step
    integer :: map_dst_xsize,map_dst_ysize
    real(r4_kind)    :: d
    logical :: found

    map_dst_xsize=size(theta_dst,1);map_dst_ysize=size(theta_dst,2)
    map_src_size =size(theta_src(:))

    do j=1,map_dst_ysize
       do i=1,map_dst_xsize
          do step = 1, map_src_size
             d = spherical_distance(theta_dst(i,j),phi_dst(i,j),theta_src(step),phi_src(step))
             if( d <= max_src_dist) then
                found = update_dest_neighbors_r4(map_src_add(i,j,:),map_src_dist(i,j,:), &
                     step,d,num_found(i,j), num_neighbors )
             endif
          enddo
       enddo
    enddo

  end subroutine full_search_r4

# 46 "horiz_interp/include/horiz_interp_spherical_r4.fh" 2
# 151 "horiz_interp/horiz_interp_spherical.F90" 2

# 1 "horiz_interp/include/horiz_interp_spherical_r8.fh" 1
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




























# 1 "horiz_interp/include/horiz_interp_spherical.inc" 1
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
  !> Initialization routine.
  !!
  !> Allocates space and initializes a derived-type variable
  !! that contains pre-computed interpolation indices and weights.
  subroutine horiz_interp_spherical_new_r8(Interp, lon_in,lat_in,lon_out,lat_out, &
       num_nbrs, max_dist, src_modulo)

    type(horiz_interp_type), intent(inout) :: Interp !< A derived type variable containing indices
                                          !! and weights for subsequent interpolations. To
                                          !! reinitialize for different grid-to-grid interpolation
                                          !! @ref horiz_interp_del must be used first.
    real(r8_kind), intent(in),  dimension(:,:)      :: lon_in !< Latitude (radians) for source data grid
    real(r8_kind), intent(in),  dimension(:,:)      :: lat_in !< Longitude (radians) for source data grid
    real(r8_kind), intent(in),  dimension(:,:)      :: lon_out !< Longitude (radians) for output data grid
    real(r8_kind), intent(in),  dimension(:,:)      :: lat_out !< Latitude (radians) for output data grid
    logical, intent(in),          optional :: src_modulo !< indicates if the boundary condition
                                          !! along zonal boundary is cyclic or not. Cyclic when true
    integer, intent(in),        optional   :: num_nbrs !< Number of nearest neighbors for regridding
                                           !! When number of neighbors within the radius max_dist
                                           !! is less than num_nbrs, All the neighbors will be used
                                           !! to interpolate onto destination grid. when number of
                                           !! neighbors within the radius max_dist is greater than
                                           !! num_nbrs, at least "num_nbrs"
  !     neighbors will be used to remap onto destination grid
    real(r8_kind), optional,             intent(in) :: max_dist !< Maximum region of influence around
                                                       !! destination grid points

    !------local variables ---------------------------------------
    integer :: i, j, n
    integer :: map_dst_xsize, map_dst_ysize, map_src_xsize, map_src_ysize
    integer :: map_src_size, num_neighbors
    real(r8_kind)    :: max_src_dist, tpi, hpi
    logical :: src_is_modulo
    real(r8_kind)    :: min_theta_dst, max_theta_dst, min_phi_dst, max_phi_dst
    real(r8_kind)    :: min_theta_src, max_theta_src, min_phi_src, max_phi_src
    integer               :: map_src_add(size(lon_out,1),size(lon_out,2),max_neighbors)
    real(r8_kind)    :: map_src_dist(size(lon_out,1),size(lon_out,2),max_neighbors)
    integer               :: num_found(size(lon_out,1),size(lon_out,2))
    integer               :: ilon(max_neighbors), jlat(max_neighbors)
    real(r8_kind), dimension(size(lon_out,1),size(lon_out,2)) :: theta_dst, phi_dst
    real(r8_kind), dimension(size(lon_in,1)*size(lon_in,2))   :: theta_src, phi_src
    integer, parameter                                             :: kindl = r8_kind

    !--------------------------------------------------------------

    pe      = mpp_pe()
    root_pe = mpp_root_pe()

    tpi = 2.0_kindl*real(PI, r8_kind); hpi = 0.5_kindl*real(PI,r8_kind)

    num_neighbors = num_nbrs_default
    if(present(num_nbrs)) num_neighbors = num_nbrs
    if (num_neighbors <= 0) call mpp_error(FATAL,'horiz_interp_spherical_mod: num_neighbors must be > 0')

    max_src_dist = real(max_dist_default, r8_kind)
    if (PRESENT(max_dist)) max_src_dist = max_dist
    Interp%horizInterpReals8_type%max_src_dist = max_src_dist

    src_is_modulo = .true.
    if (PRESENT(src_modulo)) src_is_modulo = src_modulo

    !--- check the grid size comformable
    map_dst_xsize=size(lon_out,1);map_dst_ysize=size(lon_out,2)
    map_src_xsize=size(lon_in,1); map_src_ysize=size(lon_in,2)
    map_src_size = map_src_xsize*map_src_ysize

    if (map_dst_xsize /= size(lat_out,1) .or. map_dst_ysize /= size(lat_out,2)) &
         call mpp_error(FATAL,'horiz_interp_spherical_mod: destination grids not conformable')
    if (map_src_xsize /= size(lat_in,1) .or. map_src_ysize /= size(lat_in,2)) &
         call mpp_error(FATAL,'horiz_interp_spherical_mod: source grids not conformable')

    theta_src      = reshape(lon_in,(/map_src_size/))
    phi_src        = reshape(lat_in,(/map_src_size/))
    theta_dst(:,:) = lon_out(:,:)
    phi_dst(:,:)   = lat_out(:,:)

    min_theta_dst=real(tpi, r8_kind);max_theta_dst=0.0_kindl
    min_phi_dst=real(pi, r8_kind);max_phi_dst=real(-pi, r8_kind)
    min_theta_src=real(tpi, r8_kind);max_theta_src=0.0_kindl
    min_phi_src=real(pi, r8_kind);max_phi_src=real(-pi, r8_kind)

    where(theta_dst<0.0_kindl)  theta_dst = theta_dst+real(tpi,r8_kind)

    where(theta_dst>real(tpi,r8_kind))  theta_dst = theta_dst-real(tpi,r8_kind)
    where(theta_src<0.0_kindl)  theta_src =  theta_src+real(tpi,r8_kind)
    where(theta_src>real(tpi,r8_kind))  theta_src = theta_src-real(tpi,r8_kind)

    where(phi_dst < -hpi) phi_dst = -hpi
    where(phi_dst > hpi)  phi_dst =  hpi
    where(phi_src < -hpi) phi_src = -hpi
    where(phi_src > hpi)  phi_src =  hpi

    do j=1,map_dst_ysize
       do i=1,map_dst_xsize
          min_theta_dst = min(min_theta_dst,theta_dst(i,j))
          max_theta_dst = max(max_theta_dst,theta_dst(i,j))
          min_phi_dst = min(min_phi_dst,phi_dst(i,j))
          max_phi_dst = max(max_phi_dst,phi_dst(i,j))
       enddo
    enddo

    do i=1,map_src_size
       min_theta_src = min(min_theta_src,theta_src(i))
       max_theta_src = max(max_theta_src,theta_src(i))
       min_phi_src = min(min_phi_src,phi_src(i))
       max_phi_src = max(max_phi_src,phi_src(i))
    enddo

    if (min_phi_dst < min_phi_src) print *, '=> WARNING:  latitute of dest grid exceeds src'
    if (max_phi_dst > max_phi_src) print *, '=> WARNING:  latitute of dest grid exceeds src'
    ! when src is cyclic, no need to print out the following warning.
    if(.not. src_is_modulo) then
       if (min_theta_dst < min_theta_src) print *, '=> WARNING : longitude of dest grid exceeds src'
       if (max_theta_dst > max_theta_src) print *, '=> WARNING : longitude of dest grid exceeds src'
    endif

    ! allocate memory to data type
    if(ALLOCATED(Interp%i_lon)) then
       if(size(Interp%i_lon,1) .NE. map_dst_xsize .OR.     &
          size(Interp%i_lon,2) .NE. map_dst_ysize )  call mpp_error(FATAL, &
          'horiz_interp_spherical_mod: size(Interp%i_lon(:),1) .NE. map_dst_xsize .OR. '// &
          'size(Interp%i_lon(:),2) .NE. map_dst_ysize')
    else
       allocate(Interp%i_lon(map_dst_xsize,map_dst_ysize,max_neighbors), &
            Interp%j_lat(map_dst_xsize,map_dst_ysize,max_neighbors),     &
            Interp%horizInterpReals8_type%src_dist(map_dst_xsize,map_dst_ysize,max_neighbors),  &
            Interp%num_found(map_dst_xsize,map_dst_ysize)       )
    endif

    map_src_add         = 0
    map_src_dist        = real(large, r8_kind)
    num_found           = 0

    !using radial_search to find the nearest points and corresponding distance.

    select case(trim(search_method))
    case ("radial_search") ! will be efficient, but may be not so accurate for some cases
       call radial_search(theta_src, phi_src, theta_dst, phi_dst, map_src_xsize, map_src_ysize, &
            map_src_add, map_src_dist, num_found, num_neighbors,max_src_dist,src_is_modulo)
    case ("full_search")   ! always accurate, but less efficient.
       call full_search(theta_src, phi_src, theta_dst, phi_dst, map_src_add, map_src_dist, &
            num_found, num_neighbors,max_src_dist )
    case default
       call mpp_error(FATAL,"HORIZ_INTERP_SPHERICAL_NEW_: nml search_method = "// &
                  trim(search_method)//" is not a valid namelist option")
    end select

    do j=1,map_dst_ysize
       do i=1,map_dst_xsize
          do n=1,num_found(i,j)
             if(map_src_add(i,j,n) == 0) then
                jlat(n) = 0; ilon(n) = 0
             else
                jlat(n) = map_src_add(i,j,n)/map_src_xsize + 1
                ilon(n) = map_src_add(i,j,n) - (jlat(n)-1)*map_src_xsize
                if(ilon(n) == 0) then
                   jlat(n) = jlat(n) - 1
                   ilon(n) = map_src_xsize
                endif
             endif
          enddo
          Interp%i_lon(i,j,:)         = ilon(:)
          Interp%j_lat(i,j,:)         = jlat(:)
          Interp%num_found(i,j)       = num_found(i,j)
          Interp%horizInterpReals8_type%src_dist(i,j,:)      = map_src_dist(i,j,:)
       enddo
    enddo

    Interp%nlon_src = map_src_xsize; Interp%nlat_src = map_src_ysize
    Interp%nlon_dst = map_dst_xsize; Interp%nlat_dst = map_dst_ysize
    Interp% horizInterpReals8_type % is_allocated = .true.
    Interp% interp_method = SPHERICAL

    return

  end subroutine horiz_interp_spherical_new_r8

  !#######################################################################

  !> Subroutine for performing the horizontal interpolation between two grids.
  !! horiz_interp_spherical_new_r8 must be called before calling this routine.
  subroutine horiz_interp_spherical_r8( Interp, data_in, data_out, verbose, mask_in, mask_out, missing_value)
    type(horiz_interp_type), intent(in) :: Interp !< A derived type variable containing indices
                                          !! and weights for subsequent interpolations. Returned
                                          !! by a previous call to horiz_interp_spherical_new_r8
    real(r8_kind), intent(in),  dimension(:,:)           :: data_in !< Input data on source grid
    real(r8_kind), intent(out), dimension(:,:)           :: data_out !< Output data on destination grid
    integer, intent(in),               optional :: verbose !< 0 = no output; 1 = min,max,means; 2 = most output
    real(r8_kind), intent(in), dimension(:,:),  optional :: mask_in !< Input mask, must be the same size as
                                             !! the input data. The real(r8_kind) value of mask_in must be
                                             !! in the range (0.,1.). Set mask_in=0.0 for data points
                                             !! that should not be used or have missing data
    real(r8_kind), intent(out), dimension(:,:), optional :: mask_out !< Output mask that specifies whether data
                                                                         !! was computed.
    real(r8_kind), intent(in),                  optional :: missing_value !< Used to indicate missing data

    !--- some local variables ----------------------------------------
    real(r8_kind), dimension(Interp%nlon_dst, Interp%nlat_dst,size(Interp%horizInterpReals8_type%src_dist,3)) :: wt
    real(r8_kind), dimension(Interp%nlon_src, Interp%nlat_src) :: mask_src
    real(r8_kind), dimension(Interp%nlon_dst, Interp%nlat_dst) :: mask_dst
    integer :: nlon_in, nlat_in, nlon_out, nlat_out, num_found
    integer :: m, n, i, j, k, miss_in, miss_out, i1, i2, j1, j2, iverbose
    real(r8_kind)    :: min_in, max_in, avg_in, min_out, max_out, avg_out, sum
    integer, parameter    :: kindl = r8_kind !< compiled real kind size
    !-----------------------------------------------------------------

    iverbose = 0;  if (present(verbose)) iverbose = verbose

    nlon_in  = Interp%nlon_src; nlat_in  = Interp%nlat_src
    nlon_out = Interp%nlon_dst; nlat_out = Interp%nlat_dst

    if(size(data_in,1) .ne. nlon_in .or. size(data_in,2) .ne. nlat_in ) &
         call mpp_error(FATAL,'horiz_interp_spherical_mod: size of input array incorrect')

    if(size(data_out,1) .ne. nlon_out .or. size(data_out,2) .ne. nlat_out ) &
         call mpp_error(FATAL,'horiz_interp_spherical_mod: size of output array incorrect')

    mask_src = 1.0_kindl; mask_dst = 1.0_kindl
    if(present(mask_in)) mask_src = mask_in

    do n=1,nlat_out
       do m=1,nlon_out
          ! neighbors are sorted nearest to farthest
          ! check nearest to see if it is a land point
          num_found = Interp%num_found(m,n)
          if(num_found == 0 ) then
             mask_dst(m,n) = 0.0_kindl
          else
             i1 = Interp%i_lon(m,n,1); j1 = Interp%j_lat(m,n,1)
             if (mask_src(i1,j1) .lt. 0.5_kindl ) then
                mask_dst(m,n) = 0.0_kindl
             endif

             if(num_found .gt. 1 ) then
                i2 = Interp%i_lon(m,n,2); j2 = Interp%j_lat(m,n,2)
                ! compare first 2 nearest neighbors -- if they are nearly
                ! equidistant then use this mask for robustness
                if(abs(Interp%horizInterpReals8_type%src_dist(m,n,2)-Interp%horizInterpReals8_type%src_dist(m,n,1)) .lt. &
                       real(epsln,r8_kind)) then
                   if((mask_src(i1,j1) .lt. 0.5_kindl ))  mask_dst(m,n) = 0.0_kindl
                endif
             endif

             sum=0.0_kindl
             do k=1, num_found
                if(mask_src(Interp%i_lon(m,n,k),Interp%j_lat(m,n,k)) .lt. 0.5_kindl ) then
                   wt(m,n,k) = 0.0_kindl
                else
                   if (Interp%horizInterpReals8_type%src_dist(m,n,k) <= real(epsln,r8_kind)) then
                      wt(m,n,k) = real(large, r8_kind)
                      sum = sum + real(large, r8_kind)
                   else if(Interp%horizInterpReals8_type%src_dist(m,n,k) <= Interp%horizInterpReals8_type%max_src_dist ) then
                      wt(m,n,k) = 1.0_kindl /Interp%horizInterpReals8_type%src_dist(m,n,k)
                      sum = sum+wt(m,n,k)
                   else
                      wt(m,n,k) = 0.0_kindl
                   endif
                endif
             enddo
             if (sum > real(epsln,r8_kind)) then
                do k = 1, num_found
                   wt(m,n,k) = wt(m,n,k)/sum
                enddo
             else
                mask_dst(m,n) = 0.0_kindl
             endif
          endif
       enddo
    enddo

    data_out = 0.0_kindl
    do n=1,nlat_out
       do m=1,nlon_out
          if(mask_dst(m,n) .gt. 0.5_kindl ) then
             do k=1, Interp%num_found(m,n)
                i = Interp%i_lon(m,n,k)
                j = Interp%j_lat(m,n,k)
                data_out(m,n) = data_out(m,n)+data_in(i,j)*wt(m,n,k)
             enddo
          else
             if(present(missing_value)) then
                data_out(m,n) = missing_value
             else
                data_out(m,n) = 0.0_kindl
             endif
          endif
       enddo
    enddo

    if(present(mask_out)) mask_out = mask_dst

    !***********************************************************************
    ! compute statistics: minimum, maximum, and mean
    !-----------------------------------------------------------------------

    if (iverbose > 0) then

       ! compute statistics of input data

       call stats (data_in, min_in, max_in, avg_in, miss_in, missing_value, mask=mask_src)

       ! compute statistics of output data
       call stats (data_out, min_out, max_out, avg_out, miss_out, missing_value, mask=mask_dst)

       !---- output statistics ----
       ! root_pe have the information of global mean, min and max
       if(pe == root_pe) then
          write (*,900)
          write (*,901)  min_in ,max_in, avg_in
          if (present(mask_in))  write (*,903)  miss_in
          write (*,902)  min_out,max_out,avg_out
          if (present(mask_out)) write (*,903)  miss_out
       endif
900    format (/,1x,10('-'),' output from horiz_interp ',10('-'))
901    format ('  input:  min=',f16.9,'  max=',f16.9,'  avg=',f22.15)
902    format (' output:  min=',f16.9,'  max=',f16.9,'  avg=',f22.15)
903    format ('          number of missing points = ',i6)

    endif

    return
  end subroutine horiz_interp_spherical_r8

  !#######################################################################
  !> This routine isn't used internally
  !! it's similar to the routine above, just gets the weights as an out variable
  subroutine horiz_interp_spherical_wght_r8( Interp, wt, verbose, mask_in, mask_out, missing_value)
    type (horiz_interp_type), intent(in)        :: Interp
    real(r8_kind), intent(out), dimension(:,:,:)         :: wt
    integer, intent(in),               optional :: verbose
    real(r8_kind), intent(in), dimension(:,:),  optional :: mask_in
    real(r8_kind), intent(inout), dimension(:,:), optional :: mask_out
    real(r8_kind), intent(in),                  optional :: missing_value

    !--- some local variables ----------------------------------------
    real(r8_kind), dimension(Interp%nlon_src, Interp%nlat_src) :: mask_src
    real(r8_kind), dimension(Interp%nlon_dst, Interp%nlat_dst) :: mask_dst
    integer :: nlon_in, nlat_in, nlon_out, nlat_out, num_found
    integer :: m, n, k, i1, i2, j1, j2, iverbose
    real(r8_kind)    :: sum
    integer, parameter    :: kindl = r8_kind
    !-----------------------------------------------------------------

    iverbose = 0;  if (present(verbose)) iverbose = verbose

    nlon_in  = Interp%nlon_src; nlat_in  = Interp%nlat_src
    nlon_out = Interp%nlon_dst; nlat_out = Interp%nlat_dst

    mask_src = 1.0_kindl ; mask_dst = 1.0_kindl
    if(present(mask_in)) mask_src = mask_in

    do n=1,nlat_out
       do m=1,nlon_out
          ! neighbors are sorted nearest to farthest
          ! check nearest to see if it is a land point
          num_found = Interp%num_found(m,n)

          if (num_found > num_nbrs_default) then
            if( iverbose .gt. 0) print *,'pe=',mpp_pe(),'num_found=',num_found
            num_found = num_nbrs_default
          end if

          if(num_found == 0 ) then
             mask_dst(m,n) = 0.0_kindl
          else
             i1 = Interp%i_lon(m,n,1); j1 = Interp%j_lat(m,n,1)
             if (mask_src(i1,j1) .lt. 0.5_kindl) then
                mask_dst(m,n) = 0.0_kindl
             endif

             if(num_found .gt. 1 ) then
                i2 = Interp%i_lon(m,n,2); j2 = Interp%j_lat(m,n,2)
                ! compare first 2 nearest neighbors -- if they are nearly
                ! equidistant then use this mask for robustness
                if(abs(Interp%horizInterpReals8_type%src_dist(m,n,2)-Interp%horizInterpReals8_type%src_dist(m,n,1)) .lt. &
                   real(epsln,r8_kind)) then
                   if((mask_src(i1,j1) .lt. 0.5_kindl ))  mask_dst(m,n) = 0.0_kindl
                endif
             endif

             sum=0.0_kindl
             do k=1, num_found
                if(mask_src(Interp%i_lon(m,n,k),Interp%j_lat(m,n,k)) .lt. 0.5_kindl ) then
                   wt(m,n,k) = 0.0_kindl
                else
                   if (Interp%horizInterpReals8_type%src_dist(m,n,k) <= real(epsln, r8_kind)) then
                      wt(m,n,k) = real(large, r8_kind)
                      sum = sum + real(large, r8_kind)
                   else if(Interp%horizInterpReals8_type%src_dist(m,n,k) <= Interp%horizInterpReals8_type%max_src_dist ) then
                      wt(m,n,k) = 1.0_kindl /Interp%horizInterpReals8_type%src_dist(m,n,k)
                      sum = sum+wt(m,n,k)
                   else
                      wt(m,n,k) = 0.0_kindl
                   endif
                endif
             enddo
             if (sum > real(epsln,r8_kind)) then
                do k = 1, num_found
                   wt(m,n,k) = wt(m,n,k)/sum
                enddo
             else
                mask_dst(m,n) = 0.0_kindl
             endif
          endif
       enddo
    enddo

    return
  end subroutine horiz_interp_spherical_wght_r8

  !#######################################################################

  subroutine radial_search_r8(theta_src,phi_src,theta_dst,phi_dst, map_src_xsize, map_src_ysize, &
       map_src_add, map_src_dist, num_found, num_neighbors,max_src_dist,src_is_modulo)
    real(r8_kind),    intent(in),    dimension(:)   :: theta_src, phi_src
    real(r8_kind),    intent(in),    dimension(:,:) :: theta_dst, phi_dst
    integer, intent(in)                    :: map_src_xsize, map_src_ysize
    integer, intent(out), dimension(:,:,:) :: map_src_add
    real(r8_kind),    intent(out), dimension(:,:,:) :: map_src_dist
    integer, intent(inout), dimension(:,:) :: num_found
    integer, intent(in)                    :: num_neighbors
    real(r8_kind),    intent(in)                    :: max_src_dist
    logical, intent(in)                    :: src_is_modulo

    !---------- local variables ----------------------------------------
    integer, parameter :: max_nbrs = 50
    integer :: i, j, jj, i0, j0, n, l,i_left, i_right
    integer :: map_dst_xsize, map_dst_ysize
    integer :: i_left1, i_left2, i_right1, i_right2
    integer :: map_src_size, step, step_size, bound, bound_start, bound_end
    logical :: continue_search, result, continue_radial_search
    real(r8_kind)    :: d, res
    !------------------------------------------------------------------
    map_dst_xsize=size(theta_dst,1);map_dst_ysize=size(theta_dst,2)
    map_src_size = map_src_xsize*map_src_ysize

    do j=1,map_dst_ysize
       do i=1,map_dst_xsize
          continue_search=.true.
          step = 1
          step_size = int( sqrt(real(map_src_size, kind=r8_kind )))
          do while (continue_search .and. step_size > 0)
             do while (step <= map_src_size .and. continue_search)
                ! count land points as nearest neighbors
                d = spherical_distance(theta_dst(i,j),phi_dst(i,j),theta_src(step),phi_src(step))
                if (d <= max_src_dist) then
                   result = update_dest_neighbors_r8(map_src_add(i,j,:),map_src_dist(i,j,:), &
                        step,d, num_found(i,j), num_neighbors )
                   if (result) then
                      n = 0
                      i0 = mod(step,map_src_xsize)

                      if (i0 == 0) i0 = map_src_xsize
                      res = real(step, r8_kind)/real(map_src_xsize, r8_kind)
                      j0 = ceiling(res)
                      continue_radial_search = .true.
                      do while (continue_radial_search)
                         continue_radial_search = .false.
                         n = n+1 ! radial counter
                         if(n > max_nbrs) exit
                         ! ************** left boundary *******************************
                         i_left = i0-n
                         if (i_left <= 0) then
                            if (src_is_modulo) then
                               i_left = map_src_xsize + i_left
                            else
                               i_left = 1
                            endif
                         endif

                         do l = 0, 2*n
                            jj = j0 - n - 1 + l
                            if( jj < 0) then
                               bound = ( 1 - jj )*map_src_xsize - i_left
                            else if ( jj >= map_src_ysize ) then
                               bound = ( 2*map_src_ysize - jj ) * map_src_xsize - i_left
                            else
                               bound = jj * map_src_xsize + i_left
                            endif

                            d = spherical_distance(theta_dst(i,j),phi_dst(i,j),theta_src(bound), &
                                                                 phi_src(bound))
                            if(d<=max_src_dist) then
                               result = update_dest_neighbors_r8(map_src_add(i,j,:),map_src_dist(i,j,:), &
                                    bound,d, num_found(i,j), num_neighbors)
                               if (result) continue_radial_search = .true.
                            endif
                         enddo

                         ! ***************************right boundary *******************************
                         i_right = i0+n
                         if (i_right > map_src_xsize) then
                            if (src_is_modulo) then
                               i_right = i_right - map_src_xsize
                            else
                               i_right = map_src_xsize
                            endif
                         endif

                         do l = 0, 2*n
                            jj = j0 - n - 1 + l
                            if( jj < 0) then
                               bound = ( 1 - jj )*map_src_xsize - i_right
                            else if ( jj >= map_src_ysize ) then
                               bound = ( 2*map_src_ysize - jj) * map_src_xsize - i_right

                            else
                               bound = jj * map_src_xsize + i_right
                            endif

                            d = spherical_distance(theta_dst(i,j),phi_dst(i,j),theta_src(bound), &
                                                                 phi_src(bound))
                            if(d<=max_src_dist) then
                               result = update_dest_neighbors_r8(map_src_add(i,j,:),map_src_dist(i,j,:), &
                                    bound,d, num_found(i,j), num_neighbors)
                               if (result) continue_radial_search = .true.
                            endif
                         enddo

                         ! ************************* bottom boundary **********************************
                         i_left2 = 0
                         if( i_left > i_right) then
                            i_left1 = 1
                            i_right1 = i_right
                            i_left2 = i_left
                            i_right2 = map_src_xsize
                         else
                            i_left1 = i_left
                            i_right1 = i_right
                         endif

                         jj = j0 - n - 1
                         if( jj < 0 ) then
                            bound_start = ( 1 - jj)*map_src_xsize - i_right1
                            bound_end   = ( 1 - jj)*map_src_xsize - i_left1
                         else
                            bound_start = jj * map_src_xsize + i_left1
                            bound_end = jj * map_src_xsize + i_right1
                         endif

                         bound = bound_start
                         do while (bound <= bound_end)
                            d = spherical_distance(theta_dst(i,j),phi_dst(i,j),theta_src(bound), &
                                                                 phi_src(bound))
                            if(d<=max_src_dist) then
                               result = update_dest_neighbors_r8(map_src_add(i,j,:),map_src_dist(i,j,:), &
                                    bound,d, num_found(i,j), num_neighbors)
                               if (result) continue_radial_search = .true.
                            endif
                            bound = bound + 1

                         enddo

                         if(i_left2 > 0 ) then
                            if( jj < 0 ) then
                               bound_start = ( 1 - jj)*map_src_xsize - i_right2
                               bound_end   = ( 1 - jj)*map_src_xsize - i_left2
                            else
                               bound_start = jj * map_src_xsize + i_left2
                               bound_end = jj * map_src_xsize + i_right2
                            endif

                            bound = bound_start
                            do while (bound <= bound_end)
                               d = spherical_distance(theta_dst(i,j),phi_dst(i,j),theta_src(bound), &
                                                                    phi_src(bound))
                               if(d<=max_src_dist) then
                                  result = update_dest_neighbors_r8(map_src_add(i,j,:),map_src_dist(i,j,:), &
                                       bound,d, num_found(i,j), num_neighbors)
                                  if (result) continue_radial_search = .true.
                               endif
                               bound = bound + 1
                            enddo
                         endif

                         ! ************************** top boundary ************************************
                         jj = j0 + n - 1
                         if( jj >= map_src_ysize) then
                            bound_start = ( 2*map_src_ysize - jj ) * map_src_xsize - i_right1
                            bound_end   = ( 2*map_src_ysize - jj ) * map_src_xsize - i_left1
                         else
                            bound_start = jj * map_src_xsize + i_left1
                            bound_end = jj * map_src_xsize + i_right1
                         endif

                         bound = bound_start
                         do while (bound <= bound_end)
                            d = spherical_distance(theta_dst(i,j),phi_dst(i,j),theta_src(bound), &
                                                                 phi_src(bound))
                            if(d<=max_src_dist) then
                               result = update_dest_neighbors_r8(map_src_add(i,j,:),map_src_dist(i,j,:), &
                                    bound,d, num_found(i,j), num_neighbors)
                               if (result) continue_radial_search = .true.
                            endif
                            bound = bound + 1
                         enddo

                         if(i_left2 > 0) then
                            if( jj >= map_src_ysize) then
                               bound_start = ( 2*map_src_ysize - jj ) * map_src_xsize - i_right2
                               bound_end   = ( 2*map_src_ysize - jj ) * map_src_xsize - i_left2
                            else
                               bound_start = jj * map_src_xsize + i_left2
                               bound_end = jj * map_src_xsize + i_right2
                            endif

                            bound = bound_start
                            do while (bound <= bound_end)
                               d = spherical_distance(theta_dst(i,j),phi_dst(i,j),theta_src(bound), &
                                                                    phi_src(bound))
                               if(d<=max_src_dist) then
                                  result = update_dest_neighbors_r8(map_src_add(i,j,:),map_src_dist(i,j,:), &
                                       bound,d, num_found(i,j), num_neighbors)
                                  if (result) continue_radial_search = .true.
                               endif
                               bound = bound + 1
                            enddo
                         endif

                      enddo
                      continue_search = .false. ! stop looking
                   endif
                endif
                step=step+step_size
             enddo ! search loop
             step = 1
             step_size = step_size/2
          enddo
       enddo
    enddo

    return

  end subroutine radial_search_r8


  !#####################################################################

  function update_dest_neighbors_r8(map_src_add, map_src_dist, src_add,d, num_found, min_nbrs)

    integer, intent(inout), dimension(:) :: map_src_add
    real(r8_kind), intent(inout),    dimension(:) :: map_src_dist
    integer, intent(in)                  :: src_add
    real(r8_kind), intent(in)                     :: d
    integer, intent(inout)               :: num_found
    integer, intent(in)                  :: min_nbrs

    logical :: update_dest_neighbors_r8, already_exist = .false.

    integer :: n,m

    update_dest_neighbors_r8 = .false.

    n = 0
    NLOOP : do while ( n .le. num_found )
       n = n + 1
       DIST_CHK : if (d .le. map_src_dist(n)) then
          do m=n,num_found
             if (src_add == map_src_add(m)) then
                already_exist = .true.
                exit NLOOP
             endif
          enddo
          if(num_found < max_neighbors) then
             num_found = num_found + 1
          else
             call mpp_error(FATAL,'UPDATE_DEST_NEIGHBORS_: '// &
                  'number of neighbor points found is greated than maxium neighbor points' )
          endif
          do m=num_found,n+1,-1
             map_src_add(m) = map_src_add(m-1)
             map_src_dist(m) = map_src_dist(m-1)
          enddo
          map_src_add(n) = src_add
          map_src_dist(n) = d
          update_dest_neighbors_r8 = .true.
          if( num_found > min_nbrs ) then
             if( map_src_dist(num_found) > map_src_dist(num_found-1) ) then
                num_found = num_found - 1
             endif
             if( map_src_dist(min_nbrs+1) > map_src_dist(min_nbrs) ) then
                num_found = min_nbrs
             endif
          endif
          exit NLOOP ! n loop
       endif DIST_CHK
    end do NLOOP
    if(already_exist) return

    if( .not. update_dest_neighbors_r8 ) then
       if( num_found < min_nbrs ) then
          num_found               = num_found + 1
          update_dest_neighbors_r8   = .true.
          map_src_add(num_found)  = src_add
          map_src_dist(num_found) = d
       endif
    endif


    return

  end function update_dest_neighbors_r8

  !########################################################################
!  function spherical_distance_r8(theta1,phi1,theta2,phi2)

!    real(r8_kind), intent(in) :: theta1, phi1, theta2, phi2
!    real(r8_kind) :: spherical_distance_r8

!    real(r8_kind) :: r1(3), r2(3), cross(3), s, dot, ang

    ! this is a simple, enough way to calculate distance on the sphere
    ! first, construct cartesian vectors r1 and r2
    ! then calculate the cross-product which is proportional to the area
    ! between the 2 vectors.  The angular distance is arcsin of the
    ! distancealong the sphere
    !
    ! theta is longitude and phi is latitude
    !


!    r1(1) = cos(theta1)*cos(phi1);r1(2)=sin(theta1)*cos(phi1);r1(3)=sin(phi1)
!    r2(1) = cos(theta2)*cos(phi2);r2(2)=sin(theta2)*cos(phi2);r2(3)=sin(phi2)

!    cross(1) = r1(2)*r2(3)-r1(3)*r2(2)
!    cross(2) = r1(3)*r2(1)-r1(1)*r2(3)
!    cross(3) = r1(1)*r2(2)-r1(2)*r2(1)

!    s = sqrt(cross(1)**2.+cross(2)**2.+cross(3)**2.)

!    s = min(s,real(1.0, r8_kind)-epsln)

!    dot = r1(1)*r2(1) + r1(2)*r2(2) + r1(3)*r2(3)

!    if (dot > 0) then
!       ang = asin(s)
!    else if (dot < 0) then
!       ang = pi + asin(s)  !?  original is pi - asin(s)
!    else
!       ang = pi/real(2., r8_kind)
!    endif

!    spherical_distance_r8 = abs(ang) ! in radians

!    return

!  end function spherical_distance_r8
  ! The great cycle distance
  function spherical_distance_r8(theta1,phi1,theta2,phi2)

    real(r8_kind), intent(in) :: theta1, phi1, theta2, phi2
    real(r8_kind) :: spherical_distance_r8, dot
    integer, parameter :: kindl = r8_kind

    if(theta1 == theta2 .and. phi1 == phi2) then
        spherical_distance_r8 = 0.0_kindl
        return
    endif

    dot = cos(phi1)*cos(phi2)*cos(theta1-theta2) + sin(phi1)*sin(phi2)
    if(dot > 1.0_kindl)  dot = 1.0_kindl
    if(dot < real(-1.0_kindl, r8_kind)) dot = -1.0_kindl
    spherical_distance_r8 = acos(dot)

    return

  end function spherical_distance_r8


  !#######################################################################

  subroutine full_search_r8(theta_src,phi_src,theta_dst,phi_dst,map_src_add, map_src_dist,num_found, &
                         num_neighbors,max_src_dist)
    real(r8_kind),    intent(in),    dimension(:)   :: theta_src, phi_src
    real(r8_kind),    intent(in),    dimension(:,:) :: theta_dst, phi_dst
    integer, intent(out), dimension(:,:,:) :: map_src_add
    real(r8_kind),    intent(out), dimension(:,:,:) :: map_src_dist
    integer, intent(out), dimension(:,:)   :: num_found
    integer, intent(in)                    :: num_neighbors
    real(r8_kind), intent(in)                       :: max_src_dist

    integer :: i,j,map_src_size, step
    integer :: map_dst_xsize,map_dst_ysize
    real(r8_kind)    :: d
    logical :: found

    map_dst_xsize=size(theta_dst,1);map_dst_ysize=size(theta_dst,2)
    map_src_size =size(theta_src(:))

    do j=1,map_dst_ysize
       do i=1,map_dst_xsize
          do step = 1, map_src_size
             d = spherical_distance(theta_dst(i,j),phi_dst(i,j),theta_src(step),phi_src(step))
             if( d <= max_src_dist) then
                found = update_dest_neighbors_r8(map_src_add(i,j,:),map_src_dist(i,j,:), &
                     step,d,num_found(i,j), num_neighbors )
             endif
          enddo
       enddo
    enddo

  end subroutine full_search_r8

# 46 "horiz_interp/include/horiz_interp_spherical_r8.fh" 2
# 152 "horiz_interp/horiz_interp_spherical.F90" 2

end module horiz_interp_spherical_mod
!> @}
! close documentation grouping
