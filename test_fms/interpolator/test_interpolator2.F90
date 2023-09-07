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
!! @brief unit tests for interpolator_mod
!! @author MiKyung Lee
!! @email gfdl.climate.model.info@noaa.gov
!! @description This program tests the various subroutines in interpolator_mod.  The code in
!! in test_interpolator.F90 shows how to use interpolator_mod.  This program contains the
!! actual testing.
!! TODO:

program test_interpolator2

  use fms2_io_mod, only: FmsNetcdfFile_t, UNLIMITED,                  &
                         register_field, register_variable_attribute, &
                         register_axis,                               &
                         write_data, close_file, open_file
  use mpp_mod,          only: mpp_error, FATAL, WARNING
  use time_manager_mod, only: time_type, set_calendar_type, time_manager_init, &
                              get_date_no_leap, get_date_julian, set_date, set_time, &
                              set_date_no_leap, set_date_julian, &
                              operator(+), operator(-), time_type_to_real
  use fms_mod,          only: fms_init
  use constants_mod,    only: PI
  use platform_mod,     only: r4_kind, r8_kind

  use interpolator_mod

  implicit none

  character(100), parameter :: ncfile='immadeup.o3.climatology.nc' !< fake climatology file.
  integer, parameter :: lkind=TEST_INTP_KIND_
  real(r8_kind), parameter :: tol=1.e-5_r8_kind !< the interpolation methods are not perfect.
                                                !! Will not get perfectly agreeing answers
  integer :: calendar_type

  !> climatology related variables and arrays (made up data)
  integer :: nlonlat       !< number of latitude and longitudinal center coordinates
  integer :: nlonlatb      !< number of latitude and longitudinal boundary coordinates
  integer :: ntime         !< number of time slices
  integer :: npfull        !< number of p levels
  integer :: nphalf        !< number of half p levels
  real(TEST_INTP_KIND_), allocatable :: lat(:)   !< climatology coordinates
  real(TEST_INTP_KIND_), allocatable :: lon(:)   !< climatology coordinates
  real(TEST_INTP_KIND_), allocatable :: latb(:)  !< climatology coordinates
  real(TEST_INTP_KIND_), allocatable :: lonb(:)  !< climatology coordinates
  real(TEST_INTP_KIND_), allocatable :: clim_time (:) !< climatology time
  real(TEST_INTP_KIND_), allocatable :: pfull(:) !< climatology p level
  real(TEST_INTP_KIND_), allocatable :: phalf(:) !< climatology p half level
  real(TEST_INTP_KIND_), allocatable :: ozone(:,:,:,:) !< climatology ozone data

  !> model related variables and arrays
  integer :: nlonlat_mod, nlonlatb_mod !< number of latitude and longitude coordinates in the model
  real(TEST_INTP_KIND_), allocatable :: lat_mod(:,:)  !< model coordinates
  real(TEST_INTP_KIND_), allocatable :: lon_mod(:,:)  !< model coordinates
  real(TEST_INTP_KIND_), allocatable :: latb_mod(:,:) !< model coordinates
  real(TEST_INTP_KIND_), allocatable :: lonb_mod(:,:) !< model coordinates

  type(interpolate_type) :: o3 !< recyclable interpolate_type

  logical :: test_daily_julian=.true., test_daily_noleap=.false., test_yearly_noleap=.false., test_yearly_julian=.false.
  integer :: nml_unit_var
  character(*), parameter :: nml_file='test_interpolator.nml'
  NAMELIST / test_interpolator_nml / test_daily_noleap, test_daily_julian, test_yearly_noleap, test_yearly_julian

  open(unit=nml_unit_var, file=nml_file)
  read(unit=nml_unit_var, nml=test_interpolator_nml)
  close(nml_unit_var)

  call fms_init
  call time_manager_init
  !> set data
  if(test_daily_noleap)  call set_write_data(nlonlat_in=10, nlonlat_mod_in=10, ntime_in=5, npfull_in=3, daily=.true., noleap=.true.)
  if(test_daily_julian)  call set_write_data(nlonlat_in=10, nlonlat_mod_in=10, ntime_in=5, npfull_in=3, daily=.true., noleap=.false.)
  if(test_yearly_noleap) call set_write_data(nlonlat_in=10, nlonlat_mod_in=10, ntime_in=5, npfull_in=3, yearly=.true., noleap=.true.)
  if(test_yearly_julian) call set_write_data(nlonlat_in=10, nlonlat_mod_in=10, ntime_in=5, npfull_in=3, yearly=.true., noleap=.false.)

  !> test interpolator_init with model JULIAN calendar
  calendar_type=2  !< JULIAN calendar
  call set_calendar_type(calendar_type)
  call run_test_set('JULIAN')
  !---------------------------------------------------------
  !> test interpolator_init with model NOLEAP calendar
  calendar_type=4  !< NOLEAP calendar
  call set_calendar_type(calendar_type)
  call run_test_set('NOLEAP')

  !----------------------------------------------------------------------------------
  !> Need to deallocate arrays because a new NetCDF File will be written out
  call deallocate_arrays()
  !> test interpolator_no_time_axis
  !! Write out new set of data that will have a time axis, but will have "0" time points
  !! because that's how interpolator_init is set up.
  call set_write_data(nlonlat_in=10, nlonlat_mod_in=10, ntime_in=0, npfull_in=3)
  call test_interpolator_init(o3)
  write(*,*) '                ===== test_intepolator_no_time_axis ======='
  call test_interpolator_no_time_axis(o3)

contains
  !===============================================!
  subroutine test_interpolator_init(clim_type)

    !> interopolator_init initializes the interpolate_type after reading in the
    !! climatology data.  The required arguments are clim_type, file_name, lonb_mod, latb_mod
    !! where lonb_mod and latb_mod contain the model longitude and latitude values on the grid.

    implicit none
    type(interpolate_type), intent(inout) :: clim_type
    integer, dimension(1) :: data_out_of_bounds

    data_out_of_bounds(1)=CONSTANT
    call interpolator_init(clim_type,trim(ncfile), lonb_mod, latb_mod, data_out_of_bounds=data_out_of_bounds)

  end subroutine test_interpolator_init
  !===============================================!
  subroutine test_interpolator(clim_type)

    !> call the variants of interpolator (4D-2d) that interpolates data at a given time-point
    !! The tests here do not test the "no_axis" interpolator routines
    !! This subroutine also tests obtain_interpolator_time_slices for the 2D case.

    implicit none

    type(interpolate_type), intent(inout) :: clim_type
    real(TEST_INTP_KIND_), dimension(nlonlat,nlonlat,npfull,1) :: interp_data !< last column, there is only one field
    real(TEST_INTP_KIND_), dimension(nlonlat,nlonlat,nphalf) :: phalf
    type(time_type) :: model_time, start_time
    type(time_type) :: yearly_model_time(5)
    integer :: itime, i, j, k, l

    phalf(:,:,1)=0.0000_lkind
    phalf(:,:,2)=0.0002_lkind
    phalf(:,:,3)=0.0004_lkind
    phalf(:,:,4)=0.0005_lkind

    if(test_yearly_noleap .or. test_yearly_julian) then
       yearly_model_time(1)=set_date(1850,2, 1,0,0,0)
       yearly_model_time(2)=set_date(1852,2,28,0,0,0)
       yearly_model_time(3)=set_date(1900,1, 2,0,0,0)
       yearly_model_time(4)=set_date(1905,1, 1,0,0,0)
       yearly_model_time(5)=set_date(1951,12,5,0,0,0)
    end if

    do itime=1, ntime

       !> only when clim_time is not an integer
       if(test_daily_noleap .or. test_daily_julian) model_time=get_complicated_time(clim_time(itime))

       !call print_date(model_time)

       !> test interpolator_4D_r4/8
       call interpolator(clim_type, model_time, phalf, interp_data, 'ozone')
       do i=1, npfull
          do j=1, nlonlat
             do k=1, nlonlat
                call check_answers(interp_data(k,j,i,1), ozone(k,j,i,itime), tol, 'test interpolator_4D')
             end do
          end do
       end do

       !> test interpolator_3_r4/8
       call interpolator(clim_type, model_time, phalf, interp_data(:,:,:,1), 'ozone')
       do i=1, npfull
          do j=1, nlonlat
             do k=1, nlonlat
                call check_answers(interp_data(k,j,i,1), ozone(k,j,i,itime), tol, 'test interpolator_3D')
             end do
          end do
       end do

       !> test interpolator_2D_r4/8
       call interpolator(clim_type, model_time, interp_data(:,:,1,1), 'ozone')
       do j=1, nlonlat_mod
          do k=1, nlonlat_mod
             call check_answers(interp_data(k,j,1,1), ozone(k,j,1,itime), tol, 'test interpolator_2D')
          end do
       end do

       !> Test obtain_interpolator_time_slices
       call obtain_interpolator_time_slices(clim_type,model_time)
       call interpolator(clim_type, model_time, interp_data(:,:,1,1), 'ozone')
       call unset_interpolator_time_flag(clim_type)
       do j=1, nlonlat_mod
          do k=1, nlonlat_mod
             call check_answers(interp_data(k,j,1,1), ozone(k,j,1,itime), tol, 'test interpolator_2D')
          end do
       end do

    end do

  end subroutine test_interpolator
  !===============================================!
  subroutine test_interpolator_end(clim_type)

    !> This subroutine tests interpolator_end

    implicit none

    type(interpolate_type) :: clim_type

    call interpolator_end(clim_type)

  end subroutine test_interpolator_end
  !===============================================!
  subroutine test_interpolator_no_time_axis(clim_type)

    !> This subroutine tests the variants (42-2D) of interpolator_no_time_axis

    implicit none

    type(interpolate_type) :: clim_type

    real(TEST_INTP_KIND_), dimension(nlonlat,nlonlat,nphalf-1,1) :: interp_data !< last column, there is only one field
    real(TEST_INTP_KIND_), dimension(nlonlat,nlonlat,nphalf) :: phalf
    integer :: i, j, k

    phalf(:,:,1)=0.0000_lkind
    phalf(:,:,2)=0.0002_lkind
    phalf(:,:,3)=0.0004_lkind
    phalf(:,:,4)=0.0005_lkind

    !> test interpolator_4D_no_time_axis_r4/8
    call interpolator(clim_type, phalf, interp_data, 'ozone')
    do i=1, nphalf-1
       do j=1, nlonlat
          do k=1, nlonlat
             call check_answers(interp_data(k,j,i,1), ozone(k,j,i,1), tol, 'test interpolator_4D_no_time_axis')
          end do
       end do
    end do

    !> test interpolator_3D_no_time_axis_r4/8
    call interpolator(clim_type, phalf, interp_data(:,:,:,1), 'ozone')
    do i=1, nphalf-1
       do j=1, nlonlat
          do k=1, nlonlat
             call check_answers(interp_data(k,j,i,1), ozone(k,j,i,1), tol, 'test interpolator_3D_no_time_axis')
          end do
       end do
    end do

    !> test interpolator_2D_no_time_axis_r4/8
    call interpolator(clim_type, interp_data(:,:,1,1), 'ozone')
    do j=1, nlonlat
       do k=1, nlonlat
          call check_answers(interp_data(k,j,1,1), ozone(k,j,1,1), tol, 'test interpolator_2D_no_time_axis')
       end do
    end do

  end subroutine test_interpolator_no_time_axis
  !===============================================!
  subroutine test_interpolate_type_eq

    !> This subroutine tests interpolaote_type_eq (assignment = operator)
    !! The success of "=" is insured by checking to see if interpolation with o3_copy succeds.

    implicit none

    type(interpolate_type) :: o3_copy

    o3_copy = o3
    call test_interpolator(o3_copy)

  end subroutine test_interpolate_type_eq
  !===============================================!
  subroutine test_query_interpolator

    !> This subroutne tests query_interpolator

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
  !===============================================!
  subroutine run_test_set(test_model_calendar)

    character(*), intent(in) :: test_model_calendar

    write(*,*) '                 ===== test_interpolator_init'//trim(test_model_calendar)//' ====='
    call test_interpolator_init(o3)

    !> test interpolator 2D-4D
    write(*,*) '                 ===== test_intepolator'//trim(test_model_calendar)//' ====='
    call test_interpolator(o3)

    !> test interpolate_type_eq
    !! This test has been commented out and will be included
    !! in the testing suite once fileobj cp is added into
    !! test_interpolate_type_eq
    !write(*,*) '===== test_interpolate_type_eq ====='
    !call test_interpolate_type_eq()

    !> test query_interpolator
    write(*,*) '                 ===== test_query_interpolator'//trim(test_model_calendar)//' ====='
    call test_query_interpolator()

    !> test interpolator end
    write(*,*) '                 ===== test_interpolator_end'//trim(test_model_calendar)//' ====='
    call test_interpolator_end(o3)


  end subroutine run_test_set
  !===============================================!
  function get_complicated_time(days_in)

    implicit none
    type(time_type) :: get_complicated_time
    real(TEST_INTP_KIND_), intent(in) :: days_in

    integer, parameter :: seconds_per_day=86400
    type(time_type) :: tmp_time, base_time
    integer :: yr, mo, dy, hr, mn, sc
    real(TEST_INTP_KIND_) :: frac_day

    !1849,1,1
    if(test_daily_noleap.or.test_yearly_noleap) base_time=set_date_no_leap(1849,1,1,0,0,0)
    if(test_daily_julian.or.test_yearly_julian) base_time=set_date_julian(1849,1,1,0,0,0)
    tmp_time = set_time(0,int(days_in)) + base_time
    frac_day = days_in - real(int(days_in),TEST_INTP_KIND_)
    if(test_daily_noleap .and. calendar_type==2) then
       !> daily file calendar = noleap, model calendar = julian
       call get_date_no_leap(tmp_time, yr, mo, dy, hr, mn, sc)
       get_complicated_time=set_date(yr, mo, dy, mn, sc)
    else if(test_daily_noleap .and. calendar_type==4) then
       !> daily file calendar = noleap, model calendar = noleap
       get_complicated_time=set_time(int(seconds_per_day*frac_day),int(days_in))+base_time
    else if(test_daily_julian .and. calendar_type==2) then
       !> daily file calendar = julian, model_calendar = julian
       get_complicated_time=set_time(int(seconds_per_day*frac_day),int(days_in))+base_time
    else if(test_daily_julian .and. calendar_type==4) then
       !> daily file calendar = julian, model calendar = noleap
       call get_date_julian(tmp_time, yr, mo, dy, hr, mn, sc)
       get_complicated_time=set_date(yr, mo, dy, hr, mn, sc)
    end if


  end function get_complicated_time
  !===============================================!
  subroutine check_answers(results, answers, tol, whoami)

    implicit none
    real(TEST_INTP_KIND_), intent(in) :: results, answers
    real(r8_kind), intent(in) :: tol
    character(*) :: whoami

    if (real(abs(results-answers),r8_kind).gt.tol) then
       write(*,*) '      EXPECTED ', answers, ' but computed ', results
       call mpp_error(FATAL, trim(whoami))
    end if

  end subroutine check_answers
  !===============================================!
#include "test_interpolator_write_climatology.inc"

end program test_interpolator2
