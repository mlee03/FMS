program test_interpolator

  use fms2_io_mod, only: FmsNetcdfFile_t, UNLIMITED,               &
       register_field, register_variable_attribute, register_axis, &
       write_data, close_file, open_file
  use mpp_mod, only: mpp_npes, mpp_get_current_pelist, mpp_error, FATAL
  use time_manager_mod, only: time_type, set_date, set_calendar_type, time_manager_init
  use fms_mod, only: fms_init

  use interpolator_mod

  implicit none

  character(100), parameter :: ncfile='immadeup.o3.climatology.nc'

  integer, parameter :: nlonlat=10
  integer, parameter :: nlonlatb=nlonlat+1
  integer, parameter :: calendar_type=2 !< JULIAN

  real, dimension(nlonlat,nlonlat) :: lat_mod, lon_mod
  real, dimension(nlonlat+1,nlonlat+1) :: latb_mod, lonb_mod
  type(interpolate_type) :: o3

  integer :: i,j,k

  call fms_init
  call time_manager_init
  call set_calendar_type(calendar_type)

  call write_climatology_file
  call init_model_latlon

  call test_interpolator_init
  call test_subroutine_interpolator
  call test_query_interpolator

contains
  !-----------------------------------------------!
  !-----------------------------------------------!
  !> create model lat(b) and lon(b) arrays
  subroutine init_model_latlon

    implicit none
    real :: dxy

    dxy = 1.0
    do i=1, nlonlat
       lon_mod(:,i) = real(i)*dxy
       lat_mod(i,:) = real(i)*dxy
    end do

    !!!!!!FIXXXXXXXX!!!!!
    lonb_mod(1:nlonlat,1)=0.0
    latb_mod(1,1:nlonlat)=0.0
    do i=2,nlonlat
       lonb_mod(1:nlonlat,i)=0.5*(lon_mod(:,i-1)+lon_mod(:,i))
       latb_mod(i,1:nlonlat)=0.5*(lat_mod(i-1,:)+lat_mod(i,:))
    enddo
    lonb_mod(1:nlonlat,nlonlat+1)=real(nlonlat-1)
    latb_mod(nlonlat+1,1:nlonlat)=real(nlonlat-1)

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
  subroutine test_subroutine_interpolator

    implicit none
    real, dimension(nlonlat,nlonlat) :: interp_data
    type(time_type) :: model_time

    model_time=set_date(1849,1,2)
    call interpolator(o3, model_time, interp_data, 'ozone')

  end subroutine test_subroutine_interpolator
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
