program test_interpolator

  use fms2_io_mod, only: FmsNetcdfFile_t, UNLIMITED,               &
       register_field, register_variable_attribute, register_axis, &
       write_data, close_file, open_file
  use mpp_mod, only: mpp_npes, mpp_get_current_pelist, mpp_error, FATAL
  use time_manager_mod, only: time_type, set_date, set_calendar_type, time_manager_init
  use fms_mod, only: fms_init
  use constants_mod, only: PI

  use interpolator_mod

  implicit none

  character(100), parameter :: ncfile='immadeup.o3.climatology.nc'

  integer, parameter :: nlonlat=10
  integer, parameter :: nlonlatb=nlonlat+1
  integer, parameter :: calendar_type=2 !< JULIAN
  integer, parameter :: ntime=5
  integer, parameter :: npfull=3, nphalf=npfull+1

  real, dimension(nlonlat,nlonlat) :: lat_mod, lon_mod
  real, dimension(nlonlat+1,nlonlat+1) :: latb_mod, lonb_mod
  type(interpolate_type) :: o3

  real, dimension(ntime)  :: timedata
  real, dimension(nlonlat)  :: latdata, londata
  real, dimension(nlonlatb) :: latbdata, lonbdata
  real, dimension(npfull) :: pfulldata
  real, dimension(nphalf) :: phalfdata
  real, dimension(nlonlat,nlonlat,npfull,ntime) :: ozonedata

  integer :: i,j,k

  call fms_init
  call time_manager_init
  call set_calendar_type(calendar_type)

  call write_climatology_file
  call init_model_latlon

  call test_interpolator_init
  call test_subroutine_interpolator_3D
  call test_query_interpolator

contains
  !-----------------------------------------------!
  !-----------------------------------------------!
  !> create model lat(b) and lon(b) arrays
  subroutine init_model_latlon

    implicit none
    real :: dxy

    !> the model grid and the climatology grid are the same
    do i=1, nlonlat
       lon_mod(:,i) = real(2*i-1)
       lat_mod(i,:) = real(2*i-1)
    end do

    lonb_mod(1:nlonlat,1)=0.0
    latb_mod(1,1:nlonlat)=0.0
    do i=2,nlonlat
       lonb_mod(1:nlonlat,i)=0.5*(lon_mod(:,i-1)+lon_mod(:,i))
       latb_mod(i,1:nlonlat)=0.5*(lat_mod(i-1,:)+lat_mod(i,:))
    enddo
    lonb_mod(1:nlonlat,nlonlat+1)=real(2*nlonlat-1)
    latb_mod(nlonlat+1,1:nlonlat)=real(2*nlonlat-1)

    !> convert from degrees to radians
    lon_mod=lon_mod*PI/180.0
    lat_mod=lat_mod*PI/180.0
    lonb_mod=lonb_mod*PI/180.0
    latb_mod=latb_mod*PI/180.0

  end subroutine init_model_latlon
!-----------------------------------------------!
!-----------------------------------------------!
  !> test interpolator init
  subroutine test_interpolator_init
    !> interopolator_init initializes the interpolate_type after reading in the
    !! climatology data.  The required arguments are clim_type, file_name, lonb_mod, latb_mod
    !! where lonb_mod and latb_mod contain the model longitude and latitude values on the grid.

    implicit none
    integer, dimension(1) :: data_out_of_bounds

    data_out_of_bounds(1)=CONSTANT
    call interpolator_init(o3,trim(ncfile), lonb_mod, latb_mod, data_out_of_bounds=data_out_of_bounds)

  end subroutine test_interpolator_init
!-----------------------------------------------!
!-----------------------------------------------!
  subroutine test_subroutine_interpolator_3D

    implicit none
    real, dimension(nlonlat,nlonlat,nphalf-1) :: interp_data
    real, dimension(nlonlat,nlonlat,nphalf) :: phalf
    type(time_type) :: model_time

    integer :: itime

    do itime=1, ntime

       model_time=set_date(1849,1,1+int(timedata(itime)))
       phalf(:,:,1)=0.0000
       phalf(:,:,2)=0.0002
       phalf(:,:,3)=0.0004
       phalf(:,:,4)=0.0005
       call interpolator(o3, model_time, phalf, interp_data, 'ozone')

       do i=1, nphalf-1
          do j=1, nlonlat
             do k=1, nlonlat
                if( interp_data(k,j,i).ne.ozonedata(k,j,i,itime) ) then
                   write(*,*) k,j,i,itime, interp_data(k,j,i)-ozonedata(k,j,i,itime)
                   !call mpp_error(FATAL, 'wrong')
                end if
             end do
          end do
       end do

    end do

  end subroutine test_subroutine_interpolator_3D
!-----------------------------------------------!
!-----------------------------------------------!
  subroutine test_query_interpolator

    implicit none

    character(100), parameter :: answer_field_name='ozone'
    integer, parameter :: answer_nfields=1

    character(100), allocatable :: field_names(:)
    integer :: nfields

    call query_interpolator(o3,nfields)
    if( nfields .ne. answer_nfields) call mpp_error(FATAL, '')

    allocate(field_names(nfields))
    call query_interpolator(o3,nfields,field_names)
    if(trim(answer_field_name).ne.trim(field_names(1))) call mpp_error(FATAL,'')

    deallocate(field_names)

  end subroutine test_query_interpolator
!-----------------------------------------------!
!-----------------------------------------------!
#include "test_interpolator_write_climatology.inc"

end program test_interpolator
