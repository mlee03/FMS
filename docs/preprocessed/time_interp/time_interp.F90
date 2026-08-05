# 1 "time_interp/time_interp.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "time_interp/time_interp.F90"
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
!> @defgroup time_interp_mod time_interp_mod
!! @ingroup time_interp
!! @{
!! @brief Computes a weight and dates/indices for linearly interpolating between two dates.
!! @author Bruce Wyman
!!
!! A time type is converted into two consecutive dates plus
!! a fraction representing the distance between the dates.
!! This information can be used to interpolate between the dates.
!! The dates may be expressed as years, months, or days or
!! as indices in an array.

module time_interp_mod

use time_manager_mod, only: time_type, get_date, set_date, set_time, &
                            days_in_year, days_in_month, leap_year,  &
                            time_type_to_real, real_to_time_type,    &
                            get_calendar_type, JULIAN, GREGORIAN, NO_CALENDAR, &
                            operator(+), operator(-), operator(>),   &
                            operator(<), operator( // ), operator( / ),  &
                            operator(>=), operator(<=), operator( * ), &
                            operator(==), print_date, print_time,&
                            time_list_error, date_to_string

use          fms_mod, only: write_version_number, &
                            error_mesg, FATAL, stdout, stdlog, &
                            check_nml_error, &
                            fms_error_handler
use          mpp_mod, only: input_nml_file
use     platform_mod

implicit none
private

!-----------------------------------------------------------------------

public :: time_interp_init, time_interp, fraction_of_year

!> Returns a weight and dates or indices for interpolating between two dates. The
!! interface fraction_of_year is provided for backward compatibility with the
!! previous version.\n
!!
!! Returns weight by interpolating Time between Time1 and Time2.
!! i.e. weight = (Time-Time1)/(Time2-Time1)
!! Time1 and Time2 may be specified by any of several different ways,
!! which is the reason for multiple interfaces.\n
!!
!! - If Time1 and Time2 are the begining and end of the year in which
!!   Time falls, use first interface.\n
!!
!! - If Time1 and Time2 fall on year boundaries, use second interface.\n
!!
!! - If Time1 and Time2 fall on month boundaries, use third.\n
!!
!! - If Time1 and Time2 fall on day boundaries, use fourth.\n
!!
!! - If Time1 and Time2 are consecutive elements of an assending list, use fifth.
!!   The fifth also returns the indices of Timelist between which Time falls.\n
!!
!! - The sixth interface is for cyclical data. Time_beg and Time_end specify the
!!   begining and end of a repeating period. In this case:<br>
!!              weight = (Time_adjusted - Time1) / (Time2 - Time1)
!! <br>Where:
!! @code{.F90}
!!              Time1 = Timelist(index1)
!!              Time2 = Timelist(index2)
!!              Time_adjusted = Time - N*Period
!!              Period = Time_end-Time_beg
!! @endcode
!! N is between (Time-Time_end)/Period and (Time-Time_beg)/Period
!! That is, N is the integer that results in Time_adjusted that is between Time_beg and Time_end.
!!
!! <br>Example usages:
!! @code{.F90}
!!              call time_interp( Time, weight )
!!              call time_interp( Time, weight, year1, year2 )
!!              call time_interp( Time, weight, year1, year2, month1, month2 )
!!              call time_interp( Time, weight, year1, year2, month1, month2, day1, day2 )
!!              call time_interp( Time, Timelist, weight, index1, index2 [, modtime] )
!!              call time_interp( Time, Time_beg, Time_end, Timelist, weight, index1, index2
!!              [,correct_leap_year_inconsistency])
!! @endcode
!!
!!   For all routines in this module the calendar type in module
!!   time_manager must be set.
!!
!!     The following private parameters are set by this module:
!!
!!          seconds per minute = 60
!!          minutes per hour   = 60
!!          hours   per day    = 24
!!          months  per year   = 12
!!
!! @param Time The time at which the weight is computed.
!! @param Time_beg For cyclical interpolation: Time_beg specifies the begining time of a cycle.
!! @param Time_end For cyclical interpolation: Time_end specifies the ending time of a cycle.
!! @param Timelist For cyclical interpolation: Timelist is an array of times between Time_beg and Time_end.
!!                 Must be monotonically increasing.
!! @param index1 Timelist(index1) = The largest value of Timelist which is less than mod(Time,Time_end-Time_beg)
!! @param index2 Timelist(index2) = The smallest value of Timelist which is greater than mod(Time,Time_end-Time_beg)
!! @param correct_leap_year_inconsistency Turns on a kluge for an inconsistency which may occur in a special case.
!!       When the modulo time period (i.e. Time_end - Time_beg) is a whole number of years
!!       and is not a multiple of 4, and the calendar in use has leap years, then it is
!!       likely that the interpolation will involve mapping a common year onto a leap year.
!!       In this case it is often desirable, but not absolutely necessary, to use data for
!!       Feb 28 of the leap year when it is mapped onto a common year.
!!       To turn this on, set correct_leap_year_inconsistency=.true.
!! @param weight weight = (mod(Time,Time_end-Time_beg) - Timelist(index1)) / (Timelist(index2) - Timelist(index1))
interface time_interp
    module procedure time_interp_frac_r8,  time_interp_year_r8, &
                     time_interp_month_r8, time_interp_day_r8,  &
                     time_interp_list_r8,  time_interp_modulo_r8
    module procedure time_interp_frac_r4,  time_interp_year_r4, &
                     time_interp_month_r4, time_interp_day_r4,  &
                     time_interp_list_r4,  time_interp_modulo_r4
end interface

integer, public, parameter :: NONE=0, YEAR=1, MONTH=2, DAY=3

!-----------------------------------------------------------------------

   integer, parameter ::  secmin = 60, minhour = 60, hourday = 24,  &
                          sechour = secmin*minhour,                  &
                          secday = secmin*minhour*hourday

   integer, parameter :: monyear = 12
   integer, parameter :: halfday = secday/2

   integer :: yrmod, momod, dymod
   logical :: mod_leapyear

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
# 151 "time_interp/time_interp.F90" 2

   logical :: module_is_initialized=.FALSE.
   logical :: perthlike_behavior=.FALSE.

   namelist / time_interp_nml / perthlike_behavior

contains


 subroutine time_interp_init()
   integer :: ierr, io, logunit

   if ( module_is_initialized ) return

   read (input_nml_file, time_interp_nml, iostat=io)
   ierr = check_nml_error (io, 'time_interp_nml')

   call write_version_number("TIME_INTERP_MOD", version)
   logunit = stdlog()
   write(logunit,time_interp_nml)

   module_is_initialized = .TRUE.

 end subroutine time_interp_init


!> @brief Wrapper function to return the fractional time into the current year
!! Always returns an r8_kind, conversion to r4 will be done implicitly if needed
!> @param Time time to calculate fraction with
!> @return real(kind=8) fraction of time passed in current year
 function fraction_of_year (Time)
 type(time_type), intent(in)  :: Time
 real(r8_kind) :: fraction_of_year

  call time_interp ( Time, fraction_of_year )

 end function fraction_of_year

!#######################################################################
!> Given an array of times in ascending order and a specific time returns
!! values of index1 and index2 such that the Timelist(index1)<=Time and
!! Time<=Timelist(index2), and index2=index1+1
!! index1=0, index2=1 or index=n, index2=n+1 are returned to indicate that
!! the time is out of range
subroutine bisect(Timelist,Time,index1,index2)
  type(time_type)  , intent(in)  :: Timelist(:)
  type(time_type)  , intent(in)  :: Time
  integer, optional, intent(out) :: index1, index2

  integer :: i,il,iu,n,i1,i2

  n = size(Timelist(:))

  if (Time==Timelist(1)) then
     i1 = 1 ; i2 = 2
  else if (Time==Timelist(n)) then
     i1 = n ; i2 = n+1
  else
     il = 0; iu=n+1
     do while(iu-il > 1)
        i = (iu+il)/2
        if(Timelist(i) > Time) then
           iu = i
        else
           il = i
        endif
     enddo
     i1 = il ; i2 = il+1
  endif

  if(PRESENT(index1)) index1 = i1
  if(PRESENT(index2)) index2 = i2
end subroutine bisect

!#######################################################################
!  private routines
!#######################################################################

 function year_midpt (yr)

   integer, intent(in) :: yr
   type (time_type)    :: year_midpt, year_beg, year_end


   year_beg = set_date(yr  , 1, 1)
   year_end = set_date(yr+1, 1, 1)

   year_midpt = (year_beg + year_end) / 2

 end function year_midpt

 function month_midpt (yr, mo)

   integer, intent(in) :: yr, mo
   type (time_type)    :: month_midpt, month_beg, month_end

!  --- beginning of this month ---
   month_beg = set_date(yr, mo, 1)

!  --- start of next month ---
   if (mo < 12) then
      month_end = set_date(yr, mo+1, 1)
   else
      month_end = set_date(yr+1, 1, 1)
   endif

   month_midpt = (month_beg + month_end) / 2

 end function month_midpt

function set_modtime (Tin, modtime) result (Tout)
type(time_type), intent(in) :: Tin
integer, intent(in), optional :: modtime
type(time_type)             :: Tout
integer :: yr, mo, dy, hr, mn, se, mtime

  if(present(modtime)) then
    mtime = modtime
  else
    mtime = NONE
  endif

  select case (mtime)
    case (NONE)
       Tout = Tin
    case (YEAR)
       call get_date (Tin, yr, mo, dy, hr, mn, se)
       yr = yrmod
        ! correct leap year dates
          if (.not.mod_leapyear .and. mo == 2 .and. dy > 28) then
             mo = 3; dy = dy-28
          endif
       Tout = set_date (yr, mo, dy, hr, mn, se)
    case (MONTH)
       call get_date (Tin, yr, mo, dy, hr, mn, se)
       yr = yrmod; mo = momod
       Tout = set_date (yr, mo, dy, hr, mn, se)
    case (DAY)
       call get_date (Tin, yr, mo, dy, hr, mn, se)
       yr = yrmod; mo = momod; dy = dymod
       Tout = set_date (yr, mo, dy, hr, mn, se)
  end select

end function set_modtime

subroutine error_handler(string)
  character(len=*), intent(in) :: string

  call error_mesg ('time_interp_mod', trim(string), FATAL)

end subroutine error_handler



# 1 "time_interp/include/time_interp_r4.fh" 1
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

























# 1 "time_interp/include/time_interp.inc" 1
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

 !> @brief Calculates the fractional time into the current year
 subroutine time_interp_frac_r4 ( Time, weight )

   type(time_type),   intent(in)  :: Time
   real(r4_kind), intent(out) :: weight !< fractional time

   integer         :: yr, mo, dy, hour, minute, second
   type(time_type) :: Year_beg, Year_end


   if ( .not. module_is_initialized ) call time_interp_init

!  ---- compute fractional time of year -----

     call get_date (Time, yr, mo, dy, hour, minute, second)

     Year_beg = set_date(yr  , 1, 1)
     Year_end = set_date(yr+1, 1, 1)

     weight = real( (Time - Year_beg) // (Year_end - Year_beg) , kind=r4_kind)

 end subroutine time_interp_frac_r4


 !> @brief Calculates fractional time between mid points of consecutive years
 subroutine time_interp_year_r4 ( Time, weight, year1, year2 )

   type(time_type),   intent(in)  :: Time
   real(r4_kind), intent(out) :: weight !< fractional time between midpoints of year1 and year2
   integer        ,   intent(out) :: year1, year2

   integer :: yr, mo, dy, hour, minute, second
   type (time_type) :: Mid_year, Mid_year1, Mid_year2


   if ( .not. module_is_initialized ) call time_interp_init()

      call get_date (Time, yr, mo, dy, hour, minute, second)

    ! mid point of current year
      Mid_year = year_midpt(yr)

      if ( Time >= Mid_year ) then
    ! current time is after mid point of current year
           year1  = yr
           year2  = yr+1
           Mid_year2 = year_midpt(year2)
           weight = real( (Time - Mid_year) // (Mid_year2 - Mid_year) , kind=r4_kind )
      else
    ! current time is before mid point of current year
           year2  = yr
           year1  = yr-1
           Mid_year1 = year_midpt(year1)
           weight = real( (Time - Mid_year1) // (Mid_year - Mid_year1), kind=r4_kind )
      endif

 end subroutine time_interp_year_r4

 !> @brief Calculates fractional time between mid points of consecutive months
 subroutine time_interp_month_r4 ( Time, weight, year1, year2, month1, month2 )

   type(time_type), intent(in)  :: Time
   real(r4_kind)           , intent(out) :: weight
   integer        , intent(out) :: year1, year2, month1, month2

   integer :: yr, mo, dy, hour, minute, second,  &
              mid_month, cur_month, mid1, mid2

   if ( .not. module_is_initialized ) call time_interp_init()

      call get_date (Time, yr, mo, dy, hour, minute, second)

    ! mid point of current month in seconds
      mid_month = days_in_month(Time) * halfday
    ! time into current month in seconds
      cur_month = second + secmin*minute + sechour*hour + secday*(dy-1)

      if ( cur_month >= mid_month ) then
    ! current time is after mid point of current month
           year1  = yr;  month1 = mo
           year2  = yr;  month2 = mo+1
           if (month2 > monyear)  then
              year2 = year2+1;  month2 = 1
           endif
           mid1 = mid_month
           mid2 = days_in_month(set_date(year2,month2,2)) * halfday
           weight = real(cur_month - mid1, r4_kind) / real(mid1+mid2, r4_kind)
      else
    ! current time is before mid point of current month
           year2  = yr;  month2 = mo
           year1  = yr;  month1 = mo-1
           if (month1 < 1)  then
              year1 = year1-1;  month1 = monyear
           endif
           if (year1>0) then
              mid1 = days_in_month(set_date(year1,month1,2)) * halfday
           else
              ! this can happen if we are at the beginning of year 1. In this case
              ! use December 0001 to calculate the duration of December 0000.
              ! This should work for all calendars
              mid1 = days_in_month(set_date(1,month1,2)) * halfday
           endif
           mid2 = mid_month
           weight = real(cur_month + mid1, r4_kind) / real(mid1+mid2, r4_kind)
      endif

 end subroutine time_interp_month_r4

 !> @brief Calculates fractional time between mid points of consecutive days
 subroutine time_interp_day_r4 ( Time, weight, year1, year2, month1, month2, day1, day2 )

   type(time_type), intent(in)  :: Time
   real(r4_kind), intent(out) :: weight
   integer        , intent(out) :: year1, year2, month1, month2, day1, day2

   integer :: yr, mo, dy, hour, minute, second, sday

   if ( .not. module_is_initialized ) call time_interp_init()

      call get_date (Time, yr, mo, dy, hour, minute, second)

    ! time into current day in seconds
      sday = second + secmin*minute + sechour*hour

      if ( sday >= halfday ) then
    ! current time is after mid point of day
           year1 = yr;  month1 = mo;  day1 = dy
           year2 = yr;  month2 = mo;  day2 = dy + 1
           weight  = real(sday - halfday, r4_kind) / real(secday, r4_kind)

           if (day2 > days_in_month(Time)) then
               month2 = month2 + 1
               day2 = 1
               if (month2 > monyear) then
                    month2 = 1;  year2 = year2+1
               endif
           endif
      else
    ! current time is before mid point of day
           year2 = yr;  month2 = mo   ;  day2 = dy
           year1 = yr;  month1 = mo;  day1 = dy - 1
           weight  = real(sday + halfday,r4_kind) / real(secday,r4_kind)

           if (day1 < 1) then
               month1 = month1 - 1
               if (month1 < 1) then
                   month1 = monyear;  year1 = year1-1
               endif
               day1 = days_in_month(set_date(year1,month1,2))
           endif
      endif

 end subroutine time_interp_day_r4

 !> Part of the time_interp interface, calculates for cyclical data
 !! Time_beg and Time_end mark a repeating period
 !!
 !! Finds mid points and fractional weight for a time perioid
subroutine time_interp_modulo_r4(Time, Time_beg, Time_end, Timelist, weight, index1, index2, &
                              correct_leap_year_inconsistency, err_msg)
type(time_type), intent(in)  :: Time !< a specific time value
type(time_type), intent(in)  :: Time_beg !< begining of period to search with
type(time_type), intent(in)  :: Time_end !< end of period to search with
type(time_type), intent(in)  :: Timelist(:) !< ascending time values to search between
real(r4_kind)           , intent(out) :: weight
integer        , intent(out) :: index1, index2 !< indices of bounding time values within Timelist
logical, intent(in), optional :: correct_leap_year_inconsistency!< When true turns on a kluge for an
                                !! inconsistency which may occur in a special case.
                                !! When the modulo time period (i.e. Time_end - Time_beg) is a
                                !! whole number of years and is not a multiple of 4, and the calendar
                                !! in use has leap years, then it is likely that the interpolation
                                !! will involve mapping a common year onto a leap year. In this case
                                !! it is often desirable, but not absolutely necessary, to use data
                                !! for Feb 28 of the leap year when it is mapped onto a common year.
character(len=*), intent(out), optional :: err_msg

  type(time_type) :: Period, T
  integer :: is, ie,i1,i2
  integer :: ys,ms,ds,hs,mins,ss ! components of the starting date
  integer :: ye,me,de,he,mine,se ! components of the ending date
  integer :: yt,mt,dt,ht,mint,st ! components of the current date
  integer :: dt1                 ! temporary value for day
  integer :: n                   ! size of Timelist
  integer :: stdoutunit
  logical :: correct_lyr, calendar_has_leap_years, do_the_lyr_correction
  integer, parameter :: kindl = r4_kind

  if ( .not. module_is_initialized ) call time_interp_init
  if( present(err_msg) ) err_msg = ''

  stdoutunit = stdout()
  n = size(Timelist)

  if (Time_beg>=Time_end) then
     if(fms_error_handler('time_interp_modulo', &
     'end of the specified time loop interval must be later than its beginning',err_msg)) return
  endif

  calendar_has_leap_years = (get_calendar_type() == JULIAN .or. get_calendar_type() == GREGORIAN)

  Period = Time_end-Time_beg ! period of the time axis

  if(present(correct_leap_year_inconsistency)) then
    correct_lyr = correct_leap_year_inconsistency
  else
    correct_lyr = .false.
  endif

  ! bring the requested time inside the specified time period
  T = Time

  do_the_lyr_correction = .false.

  ! Determine if the leap year correction needs to be done.
  ! It never needs to be done unless 3 conditions are met:
  ! 1) We are using a calendar with leap years
  ! 2) optional argument correct_leap_year_inconsistency is present and equals .true.
  ! 3) The modulo time period is an integer number of years
  ! If all of these are true then set do_the_lyr_correction to .true.

  if(calendar_has_leap_years .and. correct_lyr) then
    call get_date(Time_beg,ys,ms,ds,hs,mins,ss)
    call get_date(Time_end,ye,me,de,he,mine,se)
    if(ms==me.and.ds==de.and.hs==he.and.mins==mine.and.ss==se) then
      ! whole number of years
      do_the_lyr_correction = .true.
    endif
  endif

  if(do_the_lyr_correction) then
     call get_date(T,yt,mt,dt,ht,mint,st)
     yt = ys+modulo(yt-ys,ye-ys)
     dt1 = dt
     ! If it is Feb 29, but we map into a common year, use Feb 28
     if(mt==2.and.dt==29.and..not.leap_year(set_date(yt,1,1))) dt1=28
     T = set_date(yt,mt,dt1,ht,mint,st)
     if (T < Time_beg) then
       ! the requested time is within the first year,
       ! but before the starting date. So we shift it to the last year.
       if(mt==2.and.dt==29.and..not.leap_year(set_date(ye,1,1))) dt=28
       T = set_date(ye,mt,dt,ht,mint,st)
     endif
  else
     do while ( T >= Time_end )
        T = T-Period
     enddo
     do while ( T < Time_beg )
        T = T+Period
     enddo
  endif

  ! find indices of the first and last records in the Timelist that are within
  ! the requested time period.
  if (Time_end<=Timelist(1).or.Time_beg>=Timelist(n)) then
     if(get_calendar_type() == NO_CALENDAR) then
       call print_time(Time_beg,    'Time_beg'    )
       call print_time(Time_end,    'Time_end'    )
       call print_time(Timelist(1), 'Timelist(1)' )
       call print_time(Timelist(n), 'Timelist(n)' )
     else
       call print_date(Time_beg,    'Time_beg'    )
       call print_date(Time_end,    'Time_end'    )
       call print_date(Timelist(1), 'Timelist(1)' )
       call print_date(Timelist(n), 'Timelist(n)' )
     endif
     write(stdoutunit,*)'where n = size(Timelist) =',n
     if(fms_error_handler('time_interp_modulo', &
     'the entire time list is outside the specified time loop interval',err_msg)) return
  endif

  call bisect(Timelist,Time_beg,index1=i1,index2=i2)
  if (i1 < 1) then
     is = 1 ! Time_beg before lower boundary
  else if (Time_beg == Timelist(i1)) then
     is = i1 ! Time_beg right on the lower boundary
  else
     is = i2 ! Time_beg inside the interval or on upper boundary
  endif
  call bisect(Timelist,Time_end,index1=i1,index2=i2)
  if (Time_end > Timelist(i1)) then
    ie = i1
  else if (Time_end == Timelist(i1)) then
    if(Time_beg == Timelist(is)) then
      ! Timelist includes time levels at both the lower and upper ends of the period.
      ! The endpoints of Timelist specify the same point in the cycle.
      ! This ambiguity is resolved by ignoring the last time level.
      ie = i1-1
    else
      ie = i1
    endif
  else
!   This should never happen because bisect does not return i1 such that Time_end < Timelist(i1)
  endif
  if (is>=ie) then
     if(get_calendar_type() == NO_CALENDAR) then
       call print_time(Time_beg,    'Time_beg   =')
       call print_time(Time_end,    'Time_end   =')
       call print_time(Timelist(1), 'Timelist(1)=')
       call print_time(Timelist(n), 'Timelist(n)=')
     else
       call print_date(Time_beg,    'Time_beg   =')
       call print_date(Time_end,    'Time_end   =')
       call print_date(Timelist(1), 'Timelist(1)=')
       call print_date(Timelist(n), 'Timelist(n)=')
     endif
     write(stdoutunit,*)'where n = size(Timelist) =',n
     write(stdoutunit,*)'is =',is,'ie =',ie
     if(fms_error_handler('time_interp_modulo', &
     'error in calculation of time list bounds within the specified time loop interval',err_msg)) return
  endif

  ! handle special cases:
  if( T>=Timelist(ie) ) then
     ! time is after the end of the portion of the time list within the requested period
     index1 = ie;   index2 = is
     weight = real((T-Timelist(ie))//(Period-(Timelist(ie)-Timelist(is))), r4_kind )
  else if (T<Timelist(is)) then
     ! time is before the beginning of the portion of the time list within the requested period
     index1 = ie;   index2 = is
     weight = 1.0_kindl - real(((Timelist(is)-T)//(Period-(Timelist(ie)-Timelist(is)))), r4_kind )
  else
     call bisect(Timelist,T,index1,index2)
     weight = real((T-Timelist(index1)) // (Timelist(index2)-Timelist(index1)), r4_kind )
  endif

end subroutine time_interp_modulo_r4

subroutine time_interp_list_r4 ( Time, Timelist, weight, index1, index2, modtime, err_msg )
type(time_type)  , intent(in)  :: Time, Timelist(:)
real(r4_kind)             , intent(out) :: weight
integer          , intent(out) :: index1, index2
integer, optional, intent(in)  :: modtime
character(len=*), intent(out), optional :: err_msg

integer :: n, hr, mn, se, mtime
type(time_type) :: T, Ts, Te, Td, Period, Time_mod
character(len=:),allocatable :: terr, tserr, teerr
integer, parameter :: kindl = r4_kind

  if ( .not. module_is_initialized ) call time_interp_init

  if( present(err_msg) ) err_msg = ''

  weight = 0.0_kindl; index1 = 0; index2 = 0
  n = size(Timelist(:))

! setup modular time axis?
  mtime = NONE
  if (present(modtime)) then
     mtime = modtime
     Time_mod = (Timelist(1)+Timelist(n))/2
     call get_date (Time_mod, yrmod, momod, dymod, hr, mn, se)
     mod_leapyear = leap_year(Time_mod)
  endif

! set period for modulo axis
  select case (mtime)
     case (NONE)
       ! do nothing
     case (YEAR)
         Period = set_time(0,days_in_year(Time_mod))
     case (MONTH)
       ! month length must be equal
         if (days_in_month(Time_mod) /= days_in_month(Time)) then
            if(fms_error_handler ('time_interp_list','modulo months must have same length',err_msg)) return
         endif
         Period = set_time(0,days_in_month(Time_mod))
     case (DAY)
         Period = set_time(0,1)
     case default
         if(fms_error_handler ('time_interp_list','invalid value for argument modtime',err_msg)) return
  end select

! If modulo time is in effect and Timelist spans a time interval exactly equal to
! the modulo period, then the endpoints of Timelist specify the same point in the cycle.
! This ambiguity is resolved by ignoring the last time level.
  if (mtime /= NONE .and. Timelist(size(Timelist))-Timelist(1) == Period) then
     n = size(Timelist) - 1
  else
     n = size(Timelist)
  endif

! starting and ending times from list
  Ts = Timelist(1)
  Te = Timelist(n)
  Td = Te-Ts
  T  = set_modtime(Time,mtime)

! Check that Timelist does not span a time interval greater than the modulo period
  if (mtime /= NONE) then
     if (Td > Period) then
        if(fms_error_handler ('time_interp_list','period of list exceeds modulo period',err_msg)) return
     endif
  endif

! time falls on start or between start and end list values
  if ( T >= Ts .and. T < Te ) then
     call bisect(Timelist(1:n),T,index1,index2)
     weight = real( (T-Timelist(index1)) // (Timelist(index2)-Timelist(index1)), r4_kind)

! time falls before starting list value
  else if ( T < Ts ) then
     if (mtime == NONE) then
        call time_list_error(T,terr)
        call time_list_error(Ts,tserr)
        call time_list_error(Te,teerr)
        if(fms_error_handler ('time_interp_list',&
           'time '//trim(terr)//' ('//date_to_string(T)//' is before range of list '//trim(tserr)//'-'//trim(teerr)//&
           '('//date_to_string(Ts)//' - '//date_to_string(Te)//')',&
           err_msg)) return
        deallocate(terr,tserr,teerr)
     endif
     Td = Te-Ts
     weight = 1.0_kindl - real(((Ts-T) // (Period-Td)), r4_kind )
     index1 = n
     index2 = 1

! time falls on ending list value
  else if ( T == Te ) then
    if(perthlike_behavior) then
       weight = 1.0_kindl
       index1 = n-1
       index2 = n
    else
       weight = 0.0_kindl
       index1 = n
       if (mtime == NONE) then
         index2 = n
       else
         index2 = 1
       endif
    endif

! time falls after ending list value
  else if ( T > Te ) then
     if (mtime == NONE) then
        call time_list_error(T,terr)
        call time_list_error(Ts,tserr)
        call time_list_error(Te,teerr)
        if(fms_error_handler ('time_interp_list',&
           'time '//trim(terr)//' ('//date_to_string(T)//' is after range of list '//trim(tserr)//'-'//trim(teerr)//&
           '('//date_to_string(Ts)//' - '//date_to_string(Te)//')',&
           err_msg)) return
        deallocate(terr,tserr,teerr)
     endif
     Td = Te-Ts
     weight = real( (T-Te) // (Period-Td), r4_kind)
     index1 = n
     index2 = 1
  endif

end subroutine time_interp_list_r4
!> }
# 43 "time_interp/include/time_interp_r4.fh" 2
# 305 "time_interp/time_interp.F90" 2

# 1 "time_interp/include/time_interp_r8.fh" 1
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

























# 1 "time_interp/include/time_interp.inc" 1
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

 !> @brief Calculates the fractional time into the current year
 subroutine time_interp_frac_r8 ( Time, weight )

   type(time_type),   intent(in)  :: Time
   real(r8_kind), intent(out) :: weight !< fractional time

   integer         :: yr, mo, dy, hour, minute, second
   type(time_type) :: Year_beg, Year_end


   if ( .not. module_is_initialized ) call time_interp_init

!  ---- compute fractional time of year -----

     call get_date (Time, yr, mo, dy, hour, minute, second)

     Year_beg = set_date(yr  , 1, 1)
     Year_end = set_date(yr+1, 1, 1)

     weight = real( (Time - Year_beg) // (Year_end - Year_beg) , kind=r8_kind)

 end subroutine time_interp_frac_r8


 !> @brief Calculates fractional time between mid points of consecutive years
 subroutine time_interp_year_r8 ( Time, weight, year1, year2 )

   type(time_type),   intent(in)  :: Time
   real(r8_kind), intent(out) :: weight !< fractional time between midpoints of year1 and year2
   integer        ,   intent(out) :: year1, year2

   integer :: yr, mo, dy, hour, minute, second
   type (time_type) :: Mid_year, Mid_year1, Mid_year2


   if ( .not. module_is_initialized ) call time_interp_init()

      call get_date (Time, yr, mo, dy, hour, minute, second)

    ! mid point of current year
      Mid_year = year_midpt(yr)

      if ( Time >= Mid_year ) then
    ! current time is after mid point of current year
           year1  = yr
           year2  = yr+1
           Mid_year2 = year_midpt(year2)
           weight = real( (Time - Mid_year) // (Mid_year2 - Mid_year) , kind=r8_kind )
      else
    ! current time is before mid point of current year
           year2  = yr
           year1  = yr-1
           Mid_year1 = year_midpt(year1)
           weight = real( (Time - Mid_year1) // (Mid_year - Mid_year1), kind=r8_kind )
      endif

 end subroutine time_interp_year_r8

 !> @brief Calculates fractional time between mid points of consecutive months
 subroutine time_interp_month_r8 ( Time, weight, year1, year2, month1, month2 )

   type(time_type), intent(in)  :: Time
   real(r8_kind)           , intent(out) :: weight
   integer        , intent(out) :: year1, year2, month1, month2

   integer :: yr, mo, dy, hour, minute, second,  &
              mid_month, cur_month, mid1, mid2

   if ( .not. module_is_initialized ) call time_interp_init()

      call get_date (Time, yr, mo, dy, hour, minute, second)

    ! mid point of current month in seconds
      mid_month = days_in_month(Time) * halfday
    ! time into current month in seconds
      cur_month = second + secmin*minute + sechour*hour + secday*(dy-1)

      if ( cur_month >= mid_month ) then
    ! current time is after mid point of current month
           year1  = yr;  month1 = mo
           year2  = yr;  month2 = mo+1
           if (month2 > monyear)  then
              year2 = year2+1;  month2 = 1
           endif
           mid1 = mid_month
           mid2 = days_in_month(set_date(year2,month2,2)) * halfday
           weight = real(cur_month - mid1, r8_kind) / real(mid1+mid2, r8_kind)
      else
    ! current time is before mid point of current month
           year2  = yr;  month2 = mo
           year1  = yr;  month1 = mo-1
           if (month1 < 1)  then
              year1 = year1-1;  month1 = monyear
           endif
           if (year1>0) then
              mid1 = days_in_month(set_date(year1,month1,2)) * halfday
           else
              ! this can happen if we are at the beginning of year 1. In this case
              ! use December 0001 to calculate the duration of December 0000.
              ! This should work for all calendars
              mid1 = days_in_month(set_date(1,month1,2)) * halfday
           endif
           mid2 = mid_month
           weight = real(cur_month + mid1, r8_kind) / real(mid1+mid2, r8_kind)
      endif

 end subroutine time_interp_month_r8

 !> @brief Calculates fractional time between mid points of consecutive days
 subroutine time_interp_day_r8 ( Time, weight, year1, year2, month1, month2, day1, day2 )

   type(time_type), intent(in)  :: Time
   real(r8_kind), intent(out) :: weight
   integer        , intent(out) :: year1, year2, month1, month2, day1, day2

   integer :: yr, mo, dy, hour, minute, second, sday

   if ( .not. module_is_initialized ) call time_interp_init()

      call get_date (Time, yr, mo, dy, hour, minute, second)

    ! time into current day in seconds
      sday = second + secmin*minute + sechour*hour

      if ( sday >= halfday ) then
    ! current time is after mid point of day
           year1 = yr;  month1 = mo;  day1 = dy
           year2 = yr;  month2 = mo;  day2 = dy + 1
           weight  = real(sday - halfday, r8_kind) / real(secday, r8_kind)

           if (day2 > days_in_month(Time)) then
               month2 = month2 + 1
               day2 = 1
               if (month2 > monyear) then
                    month2 = 1;  year2 = year2+1
               endif
           endif
      else
    ! current time is before mid point of day
           year2 = yr;  month2 = mo   ;  day2 = dy
           year1 = yr;  month1 = mo;  day1 = dy - 1
           weight  = real(sday + halfday,r8_kind) / real(secday,r8_kind)

           if (day1 < 1) then
               month1 = month1 - 1
               if (month1 < 1) then
                   month1 = monyear;  year1 = year1-1
               endif
               day1 = days_in_month(set_date(year1,month1,2))
           endif
      endif

 end subroutine time_interp_day_r8

 !> Part of the time_interp interface, calculates for cyclical data
 !! Time_beg and Time_end mark a repeating period
 !!
 !! Finds mid points and fractional weight for a time perioid
subroutine time_interp_modulo_r8(Time, Time_beg, Time_end, Timelist, weight, index1, index2, &
                              correct_leap_year_inconsistency, err_msg)
type(time_type), intent(in)  :: Time !< a specific time value
type(time_type), intent(in)  :: Time_beg !< begining of period to search with
type(time_type), intent(in)  :: Time_end !< end of period to search with
type(time_type), intent(in)  :: Timelist(:) !< ascending time values to search between
real(r8_kind)           , intent(out) :: weight
integer        , intent(out) :: index1, index2 !< indices of bounding time values within Timelist
logical, intent(in), optional :: correct_leap_year_inconsistency!< When true turns on a kluge for an
                                !! inconsistency which may occur in a special case.
                                !! When the modulo time period (i.e. Time_end - Time_beg) is a
                                !! whole number of years and is not a multiple of 4, and the calendar
                                !! in use has leap years, then it is likely that the interpolation
                                !! will involve mapping a common year onto a leap year. In this case
                                !! it is often desirable, but not absolutely necessary, to use data
                                !! for Feb 28 of the leap year when it is mapped onto a common year.
character(len=*), intent(out), optional :: err_msg

  type(time_type) :: Period, T
  integer :: is, ie,i1,i2
  integer :: ys,ms,ds,hs,mins,ss ! components of the starting date
  integer :: ye,me,de,he,mine,se ! components of the ending date
  integer :: yt,mt,dt,ht,mint,st ! components of the current date
  integer :: dt1                 ! temporary value for day
  integer :: n                   ! size of Timelist
  integer :: stdoutunit
  logical :: correct_lyr, calendar_has_leap_years, do_the_lyr_correction
  integer, parameter :: kindl = r8_kind

  if ( .not. module_is_initialized ) call time_interp_init
  if( present(err_msg) ) err_msg = ''

  stdoutunit = stdout()
  n = size(Timelist)

  if (Time_beg>=Time_end) then
     if(fms_error_handler('time_interp_modulo', &
     'end of the specified time loop interval must be later than its beginning',err_msg)) return
  endif

  calendar_has_leap_years = (get_calendar_type() == JULIAN .or. get_calendar_type() == GREGORIAN)

  Period = Time_end-Time_beg ! period of the time axis

  if(present(correct_leap_year_inconsistency)) then
    correct_lyr = correct_leap_year_inconsistency
  else
    correct_lyr = .false.
  endif

  ! bring the requested time inside the specified time period
  T = Time

  do_the_lyr_correction = .false.

  ! Determine if the leap year correction needs to be done.
  ! It never needs to be done unless 3 conditions are met:
  ! 1) We are using a calendar with leap years
  ! 2) optional argument correct_leap_year_inconsistency is present and equals .true.
  ! 3) The modulo time period is an integer number of years
  ! If all of these are true then set do_the_lyr_correction to .true.

  if(calendar_has_leap_years .and. correct_lyr) then
    call get_date(Time_beg,ys,ms,ds,hs,mins,ss)
    call get_date(Time_end,ye,me,de,he,mine,se)
    if(ms==me.and.ds==de.and.hs==he.and.mins==mine.and.ss==se) then
      ! whole number of years
      do_the_lyr_correction = .true.
    endif
  endif

  if(do_the_lyr_correction) then
     call get_date(T,yt,mt,dt,ht,mint,st)
     yt = ys+modulo(yt-ys,ye-ys)
     dt1 = dt
     ! If it is Feb 29, but we map into a common year, use Feb 28
     if(mt==2.and.dt==29.and..not.leap_year(set_date(yt,1,1))) dt1=28
     T = set_date(yt,mt,dt1,ht,mint,st)
     if (T < Time_beg) then
       ! the requested time is within the first year,
       ! but before the starting date. So we shift it to the last year.
       if(mt==2.and.dt==29.and..not.leap_year(set_date(ye,1,1))) dt=28
       T = set_date(ye,mt,dt,ht,mint,st)
     endif
  else
     do while ( T >= Time_end )
        T = T-Period
     enddo
     do while ( T < Time_beg )
        T = T+Period
     enddo
  endif

  ! find indices of the first and last records in the Timelist that are within
  ! the requested time period.
  if (Time_end<=Timelist(1).or.Time_beg>=Timelist(n)) then
     if(get_calendar_type() == NO_CALENDAR) then
       call print_time(Time_beg,    'Time_beg'    )
       call print_time(Time_end,    'Time_end'    )
       call print_time(Timelist(1), 'Timelist(1)' )
       call print_time(Timelist(n), 'Timelist(n)' )
     else
       call print_date(Time_beg,    'Time_beg'    )
       call print_date(Time_end,    'Time_end'    )
       call print_date(Timelist(1), 'Timelist(1)' )
       call print_date(Timelist(n), 'Timelist(n)' )
     endif
     write(stdoutunit,*)'where n = size(Timelist) =',n
     if(fms_error_handler('time_interp_modulo', &
     'the entire time list is outside the specified time loop interval',err_msg)) return
  endif

  call bisect(Timelist,Time_beg,index1=i1,index2=i2)
  if (i1 < 1) then
     is = 1 ! Time_beg before lower boundary
  else if (Time_beg == Timelist(i1)) then
     is = i1 ! Time_beg right on the lower boundary
  else
     is = i2 ! Time_beg inside the interval or on upper boundary
  endif
  call bisect(Timelist,Time_end,index1=i1,index2=i2)
  if (Time_end > Timelist(i1)) then
    ie = i1
  else if (Time_end == Timelist(i1)) then
    if(Time_beg == Timelist(is)) then
      ! Timelist includes time levels at both the lower and upper ends of the period.
      ! The endpoints of Timelist specify the same point in the cycle.
      ! This ambiguity is resolved by ignoring the last time level.
      ie = i1-1
    else
      ie = i1
    endif
  else
!   This should never happen because bisect does not return i1 such that Time_end < Timelist(i1)
  endif
  if (is>=ie) then
     if(get_calendar_type() == NO_CALENDAR) then
       call print_time(Time_beg,    'Time_beg   =')
       call print_time(Time_end,    'Time_end   =')
       call print_time(Timelist(1), 'Timelist(1)=')
       call print_time(Timelist(n), 'Timelist(n)=')
     else
       call print_date(Time_beg,    'Time_beg   =')
       call print_date(Time_end,    'Time_end   =')
       call print_date(Timelist(1), 'Timelist(1)=')
       call print_date(Timelist(n), 'Timelist(n)=')
     endif
     write(stdoutunit,*)'where n = size(Timelist) =',n
     write(stdoutunit,*)'is =',is,'ie =',ie
     if(fms_error_handler('time_interp_modulo', &
     'error in calculation of time list bounds within the specified time loop interval',err_msg)) return
  endif

  ! handle special cases:
  if( T>=Timelist(ie) ) then
     ! time is after the end of the portion of the time list within the requested period
     index1 = ie;   index2 = is
     weight = real((T-Timelist(ie))//(Period-(Timelist(ie)-Timelist(is))), r8_kind )
  else if (T<Timelist(is)) then
     ! time is before the beginning of the portion of the time list within the requested period
     index1 = ie;   index2 = is
     weight = 1.0_kindl - real(((Timelist(is)-T)//(Period-(Timelist(ie)-Timelist(is)))), r8_kind )
  else
     call bisect(Timelist,T,index1,index2)
     weight = real((T-Timelist(index1)) // (Timelist(index2)-Timelist(index1)), r8_kind )
  endif

end subroutine time_interp_modulo_r8

subroutine time_interp_list_r8 ( Time, Timelist, weight, index1, index2, modtime, err_msg )
type(time_type)  , intent(in)  :: Time, Timelist(:)
real(r8_kind)             , intent(out) :: weight
integer          , intent(out) :: index1, index2
integer, optional, intent(in)  :: modtime
character(len=*), intent(out), optional :: err_msg

integer :: n, hr, mn, se, mtime
type(time_type) :: T, Ts, Te, Td, Period, Time_mod
character(len=:),allocatable :: terr, tserr, teerr
integer, parameter :: kindl = r8_kind

  if ( .not. module_is_initialized ) call time_interp_init

  if( present(err_msg) ) err_msg = ''

  weight = 0.0_kindl; index1 = 0; index2 = 0
  n = size(Timelist(:))

! setup modular time axis?
  mtime = NONE
  if (present(modtime)) then
     mtime = modtime
     Time_mod = (Timelist(1)+Timelist(n))/2
     call get_date (Time_mod, yrmod, momod, dymod, hr, mn, se)
     mod_leapyear = leap_year(Time_mod)
  endif

! set period for modulo axis
  select case (mtime)
     case (NONE)
       ! do nothing
     case (YEAR)
         Period = set_time(0,days_in_year(Time_mod))
     case (MONTH)
       ! month length must be equal
         if (days_in_month(Time_mod) /= days_in_month(Time)) then
            if(fms_error_handler ('time_interp_list','modulo months must have same length',err_msg)) return
         endif
         Period = set_time(0,days_in_month(Time_mod))
     case (DAY)
         Period = set_time(0,1)
     case default
         if(fms_error_handler ('time_interp_list','invalid value for argument modtime',err_msg)) return
  end select

! If modulo time is in effect and Timelist spans a time interval exactly equal to
! the modulo period, then the endpoints of Timelist specify the same point in the cycle.
! This ambiguity is resolved by ignoring the last time level.
  if (mtime /= NONE .and. Timelist(size(Timelist))-Timelist(1) == Period) then
     n = size(Timelist) - 1
  else
     n = size(Timelist)
  endif

! starting and ending times from list
  Ts = Timelist(1)
  Te = Timelist(n)
  Td = Te-Ts
  T  = set_modtime(Time,mtime)

! Check that Timelist does not span a time interval greater than the modulo period
  if (mtime /= NONE) then
     if (Td > Period) then
        if(fms_error_handler ('time_interp_list','period of list exceeds modulo period',err_msg)) return
     endif
  endif

! time falls on start or between start and end list values
  if ( T >= Ts .and. T < Te ) then
     call bisect(Timelist(1:n),T,index1,index2)
     weight = real( (T-Timelist(index1)) // (Timelist(index2)-Timelist(index1)), r8_kind)

! time falls before starting list value
  else if ( T < Ts ) then
     if (mtime == NONE) then
        call time_list_error(T,terr)
        call time_list_error(Ts,tserr)
        call time_list_error(Te,teerr)
        if(fms_error_handler ('time_interp_list',&
           'time '//trim(terr)//' ('//date_to_string(T)//' is before range of list '//trim(tserr)//'-'//trim(teerr)//&
           '('//date_to_string(Ts)//' - '//date_to_string(Te)//')',&
           err_msg)) return
        deallocate(terr,tserr,teerr)
     endif
     Td = Te-Ts
     weight = 1.0_kindl - real(((Ts-T) // (Period-Td)), r8_kind )
     index1 = n
     index2 = 1

! time falls on ending list value
  else if ( T == Te ) then
    if(perthlike_behavior) then
       weight = 1.0_kindl
       index1 = n-1
       index2 = n
    else
       weight = 0.0_kindl
       index1 = n
       if (mtime == NONE) then
         index2 = n
       else
         index2 = 1
       endif
    endif

! time falls after ending list value
  else if ( T > Te ) then
     if (mtime == NONE) then
        call time_list_error(T,terr)
        call time_list_error(Ts,tserr)
        call time_list_error(Te,teerr)
        if(fms_error_handler ('time_interp_list',&
           'time '//trim(terr)//' ('//date_to_string(T)//' is after range of list '//trim(tserr)//'-'//trim(teerr)//&
           '('//date_to_string(Ts)//' - '//date_to_string(Te)//')',&
           err_msg)) return
        deallocate(terr,tserr,teerr)
     endif
     Td = Te-Ts
     weight = real( (T-Te) // (Period-Td), r8_kind)
     index1 = n
     index2 = 1
  endif

end subroutine time_interp_list_r8
!> }
# 43 "time_interp/include/time_interp_r8.fh" 2
# 306 "time_interp/time_interp.F90" 2

end module time_interp_mod

!> @}
! close documentation grouping
