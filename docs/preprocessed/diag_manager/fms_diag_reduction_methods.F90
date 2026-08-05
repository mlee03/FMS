# 1 "diag_manager/fms_diag_reduction_methods.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "diag_manager/fms_diag_reduction_methods.F90"
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

!> @defgroup fms_diag_reduction_methods_mod fms_diag_reduction_methods_mod
!! @ingroup diag_manager
!! @{
!! @brief fms_diag_reduction_methods_mod contains routines that are meant to be used for
!! error checking and setting up to do the reduction methods

module fms_diag_reduction_methods_mod
  use platform_mod, only: r8_kind, r4_kind
  use fms_diag_bbox_mod, only: fmsDiagIbounds_type
  use fms_string_utils_mod, only: string
  use diag_data_mod, only: time_diurnal, time_rms
  use mpp_mod
  implicit none
  private

  public :: check_indices_order, init_mask, set_weight
  public :: do_time_none, do_time_min, do_time_max, do_time_sum_update, time_update_done

  !> @brief Does the time_none reduction method. See include/fms_diag_reduction_methods.inc
  !TODO This needs to be extended to integers
  interface do_time_none
    module procedure do_time_none_r4, do_time_none_r8
  end interface do_time_none

  !> @brief Does the time_min reduction method. See include/fms_diag_reduction_methods.inc
  !TODO This needs to be extended to integers
  interface do_time_min
    module procedure do_time_min_r4, do_time_min_r8
  end interface do_time_min

  !> @brief Does the time_max reduction method. See include/fms_diag_reduction_methods.inc
  !TODO This needs to be extended to integers
  interface do_time_max
    module procedure do_time_max_r4, do_time_max_r8
  end interface do_time_max

  !> @brief Sum update updates the buffer for any reductions that involve summation
  !! (ie. time_sum, avg, rms, pow)
  !!TODO This needs to be extended to integers
  interface do_time_sum_update
    module procedure do_time_sum_update_r4, do_time_sum_update_r8
  end interface

  !> @brief Finishes a reduction that involves an average
  !! (ie. time_avg, rms, pow)
  !! This takes the average at the end of the time step
  interface time_update_done
    module procedure sum_update_done_r4, sum_update_done_r8
  end interface

  !> @brief Updates the buffer for any reductions that involve summation
  !! (ie. time_sum, avg, rms, pow)
  !! In this case the mask is present
  interface sum_mask
    module procedure sum_mask_r4, sum_mask_r8
  end interface

  !> @brief Updates the buffer for any reductions that involve summation
  !! (ie. time_sum, avg, rms, pow)
  !! In this case the mask is present and it varies over time
  interface sum_mask_variant
    module procedure sum_mask_variant_r4, sum_mask_variant_r8
  end interface sum_mask_variant

  !> @brief Updates the buffer for any reductions that involve summation
  !! (ie. time_sum, avg, rms, pow)
  !! In this case the mask is not present
  interface sum_no_mask
    module procedure sum_no_mask_r4, sum_no_mask_r8
  end interface sum_no_mask

  contains

  !> @brief Checks improper combinations of is, ie, js, and je.
  !! @return The error message, empty string if no errors were found
  !> @note accept_data works in either one or another of two modes.
  !! 1. Input field is a window (e.g. FMS physics)
  !! 2. Input field includes halo data
  !! It cannot handle a window of data that has halos.
  !! (A field with no windows or halos can be thought of as a special case of either mode.)
  !! The logic for indexing is quite different for these two modes, but is not clearly separated.
  !! If both the beggining and ending indices are present, then field is assumed to have halos.
  !! If only beggining indices are present, then field is assumed to be a window.
  !> @par
  !! There are a number of ways a user could mess up this logic, depending on the combination
  !! of presence/absence of is,ie,js,je. The checks below should catch improper combinations.
  pure function check_indices_order(is_in, ie_in, js_in, je_in) &
  result(error_msg)
    integer, intent(in), optional :: is_in, ie_in, js_in, je_in !< Indices passed to fms_diag_accept_data()
    character(len=128) :: error_msg !< An error message used only for testing purpose!!!

    error_msg = ""
    IF ( PRESENT(ie_in) ) THEN
      IF ( .NOT.PRESENT(is_in) ) THEN
        error_msg = 'ie_in present without is_in'
        return
      END IF
      IF ( PRESENT(js_in) .AND. .NOT.PRESENT(je_in) ) THEN
        error_msg = 'is_in and ie_in present, but js_in present without je_in'
        return
      END IF
    END IF

    IF ( PRESENT(je_in) ) THEN
      IF ( .NOT.PRESENT(js_in) ) THEN
        error_msg = 'je_in present without js_in'
        return
      END IF
      IF ( PRESENT(is_in) .AND. .NOT.PRESENT(ie_in) ) THEN
        error_msg = 'js_in and je_in present, but is_in present without ie_in'
        return
      END IF
    END IF
  end function check_indices_order

  !> @brief Sets the logical mask based on mask or rmask
  !> @return logical mask
  function init_mask(rmask, mask, field) &
  result(oor_mask)
    LOGICAL,  DIMENSION(:,:,:,:), allocatable, INTENT(in) :: mask  !< The location of the mask
    CLASS(*), DIMENSION(:,:,:,:), allocatable, INTENT(in) :: rmask !< The masking values
    CLASS(*), DIMENSION(:,:,:,:),          intent(in) :: field !< Field_data

    logical, allocatable, dimension(:,:,:,:) :: oor_mask !< mask

    ALLOCATE(oor_mask(SIZE(field, 1), SIZE(field, 2), SIZE(field, 3), SIZE(field, 4)))
    oor_mask = .true.

    if (allocated(mask)) then
      oor_mask = mask
    elseif (allocated(rmask)) then
      select type (rmask)
      type is (real(kind=r8_kind))
        WHERE (rmask < 0.5_r8_kind) oor_mask = .FALSE.
      type is (real(kind=r4_kind))
        WHERE (rmask < 0.5_r4_kind) oor_mask = .FALSE.
      end select
    endif

  end function init_mask

  !> @brief Sets the weight based on the weight passed into send_data (1.0_r8_kind if the weight is not passed in)
  !! The weight will be saved as an r8 and converted to r4 as needed
  !! @return weight to use when averaging
  pure function set_weight(weight) &
  result(out_weight)
    CLASS(*), INTENT(in), OPTIONAL :: weight !< The weight use when averaging

    real(kind=r8_kind) :: out_weight

    out_weight = 1.0_r8_kind
    if (present(weight)) then
      select type(weight)
      type is (real(kind=r8_kind))
        out_weight = real(weight, kind = r8_kind)
      type is (real(kind=r4_kind))
        out_Weight = real(weight, kind = r8_kind)
      end select
    endif
  end function set_weight


# 1 "diag_manager/include/fms_diag_reduction_methods_r4.fh" 1
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





























# 1 "diag_manager/include/fms_diag_reduction_methods.inc" 1
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

! for any debug prints




!> @brief Do the time_none reduction method (i.e copy the correct portion of the input data)
subroutine do_time_none_r4 (data_out, data_in, mask, is_masked, bounds_in, bounds_out, missing_value)
  real(r4_kind),       intent(inout) :: data_out(:,:,:,:,:) !< output data
  real(r4_kind),       intent(in)    :: data_in(:,:,:,:)    !< data to update the buffer with
  logical,                   intent(in)    :: mask(:,:,:,:)       !< mask
  logical,                   intent(in)    :: is_masked           !< .True. if the field is using a mask
  type(fmsDiagIbounds_type), intent(in)    :: bounds_in           !< indices indicating the correct portion
                                                                  !! of the input buffer
  type(fmsDiagIbounds_type), intent(in)    :: bounds_out          !< indices indicating the correct portion
                                                                  !! of the output buffer
  real(r4_kind),       intent(in)    :: missing_value       !< Missing_value for data points that are masked

  integer :: is_in, ie_in, js_in, je_in, ks_in, ke_in       !< Starting and ending indices of each dimention for
                                                            !! the input buffer
  integer :: is_out, ie_out, js_out, je_out, ks_out, ke_out !< Starting and ending indices of each dimention for
                                                            !! the output buffer

  is_out = bounds_out%get_imin()
  ie_out = bounds_out%get_imax()
  js_out = bounds_out%get_jmin()
  je_out = bounds_out%get_jmax()
  ks_out = bounds_out%get_kmin()
  ke_out = bounds_out%get_kmax()

  is_in = bounds_in%get_imin()
  ie_in = bounds_in%get_imax()
  js_in = bounds_in%get_jmin()
  je_in = bounds_in%get_jmax()
  ks_in = bounds_in%get_kmin()
  ke_in = bounds_in%get_kmax()

  if (is_masked) then
    where (mask(is_in:ie_in, js_in:je_in, ks_in:ke_in, :))
      data_out(is_out:ie_out, js_out:je_out, ks_out:ke_out, :, 1) = &
      data_in(is_in:ie_in, js_in:je_in, ks_in:ke_in, :)
    elsewhere
      data_out(is_out:ie_out, js_out:je_out, ks_out:ke_out, :, 1) = missing_value
    end where
  else
    data_out(is_out:ie_out, js_out:je_out, ks_out:ke_out, :, 1) = &
      data_in(is_in:ie_in, js_in:je_in, ks_in:ke_in, :)
  endif

end subroutine do_time_none_r4

!> @brief Do the time_min reduction method (i.e maintain the minimum value of the averaging time)
subroutine do_time_min_r4 (data_out, data_in, mask, is_masked, bounds_in, bounds_out, missing_value)
  real(r4_kind),       intent(inout) :: data_out(:,:,:,:,:) !< output data
  real(r4_kind),       intent(in)    :: data_in(:,:,:,:)    !< data to update the buffer with
  logical,                   intent(in)    :: mask(:,:,:,:)       !< mask
  logical,                   intent(in)    :: is_masked           !< .True. if the field is using a mask
  type(fmsDiagIbounds_type), intent(in)    :: bounds_in           !< indices indicating the correct portion
                                                                  !! of the input buffer
  type(fmsDiagIbounds_type), intent(in)    :: bounds_out          !< indices indicating the correct portion
                                                                  !! of the output buffer
  real(r4_kind),       intent(in)    :: missing_value       !< Missing_value for data points that are masked

  integer :: is_in, ie_in, js_in, je_in, ks_in, ke_in       !< Starting and ending indices of each dimention for
                                                            !! the input buffer
  integer :: is_out, ie_out, js_out, je_out, ks_out, ke_out !< Starting and ending indices of each dimention for
                                                            !! the output buffer

  integer :: i, j, k, l !< For looping

  is_out = bounds_out%get_imin()
  ie_out = bounds_out%get_imax()
  js_out = bounds_out%get_jmin()
  je_out = bounds_out%get_jmax()
  ks_out = bounds_out%get_kmin()
  ke_out = bounds_out%get_kmax()

  is_in = bounds_in%get_imin()
  ie_in = bounds_in%get_imax()
  js_in = bounds_in%get_jmin()
  je_in = bounds_in%get_jmax()
  ks_in = bounds_in%get_kmin()
  ke_in = bounds_in%get_kmax()

  !> Separated this loops for performance. If is_masked = .false. (i.e "mask" and "rmask" were never passed in)
  !! then mask will always be .True. so the if (mask) is redudant.
  if (is_masked) then
    do l = 0, size(data_out, 4) - 1
      do k = 0, ke_out - ks_out
        do j = 0, je_out - js_out
          do i = 0, ie_out - is_out
            if (mask(is_in + i, js_in + j, ks_in + k, l + 1)) then
              if (data_out(is_out + i, js_out + j, ks_out + k, l + 1, 1) .gt. &
                data_in(is_in + i, js_in + j, ks_in + k, l + 1) ) then
                  data_out(is_out + i, js_out + j, ks_out + k, l + 1, 1) = &
                    data_in(is_in +i, js_in + j, ks_in + k, l + 1)
              endif
            else
              data_out(is_out + i, js_out + j, ks_out + k, l + 1, 1) = missing_value
            endif
          enddo
        enddo
      enddo
    enddo
  else
    do l = 0, size(data_out, 4) - 1
      do k = 0, ke_out - ks_out
        do j = 0, je_out - js_out
          do i = 0, ie_out - is_out
            if (data_out(is_out + i, js_out + j, ks_out + k, l + 1, 1) .gt. &
              data_in(is_in + i, js_in + j, ks_in + k, l + 1) ) then
                data_out(is_out + i, js_out + j, ks_out + k, l + 1, 1) = &
                  data_in(is_in +i, js_in + j, ks_in + k, l + 1)
            endif
          enddo
        enddo
      enddo
    enddo
  endif

end subroutine do_time_min_r4

!> @brief Do the time_max reduction method (i.e maintain the maximum value of the averaging time)
subroutine do_time_max_r4 (data_out, data_in, mask, is_masked, bounds_in, bounds_out, missing_value)
  real(r4_kind),       intent(inout) :: data_out(:,:,:,:,:) !< output data
  real(r4_kind),       intent(in)    :: data_in(:,:,:,:)    !< data to update the buffer with
  logical,                   intent(in)    :: mask(:,:,:,:)       !< mask
  logical,                   intent(in)    :: is_masked           !< .True. if the field is using a mask
  type(fmsDiagIbounds_type), intent(in)    :: bounds_in           !< indices indicating the correct portion
                                                                  !! of the input buffer
  type(fmsDiagIbounds_type), intent(in)    :: bounds_out          !< indices indicating the correct portion
                                                                  !! of the output buffer
  real(r4_kind),       intent(in)    :: missing_value       !< Missing_value for data points that are masked

  integer :: is_in, ie_in, js_in, je_in, ks_in, ke_in       !< Starting and ending indices of each dimention for
                                                            !! the input buffer
  integer :: is_out, ie_out, js_out, je_out, ks_out, ke_out !< Starting and ending indices of each dimention for
                                                            !! the output buffer

  integer :: i, j, k, l !< For looping

  is_out = bounds_out%get_imin()
  ie_out = bounds_out%get_imax()
  js_out = bounds_out%get_jmin()
  je_out = bounds_out%get_jmax()
  ks_out = bounds_out%get_kmin()
  ke_out = bounds_out%get_kmax()

  is_in = bounds_in%get_imin()
  ie_in = bounds_in%get_imax()
  js_in = bounds_in%get_jmin()
  je_in = bounds_in%get_jmax()
  ks_in = bounds_in%get_kmin()
  ke_in = bounds_in%get_kmax()

  !> Separated this loops for performance. If is_masked = .false. (i.e "mask" and "rmask" were never passed in)
  !! then mask will always be .True. so the if (mask) is redudant.
  if (is_masked) then
    do l = 0, size(data_out, 4) - 1
      do k = 0, ke_out - ks_out
        do j = 0, je_out - js_out
          do i = 0, ie_out - is_out
            if (mask(is_in + i, js_in + j, ks_in + k, l + 1)) then
              if (data_out(is_out + i, js_out + j, ks_out + k, l + 1, 1) .lt. &
                data_in(is_in + i, js_in + j, ks_in + k, l + 1) ) then
                  data_out(is_out + i, js_out + j, ks_out + k, l + 1, 1) = &
                    data_in(is_in +i, js_in + j, ks_in + k, l + 1)
              endif
            else
              data_out(is_out + i, js_out + j, ks_out + k, l + 1, 1) = missing_value
            endif
          enddo
        enddo
      enddo
    enddo
  else
    do l = 0, size(data_out, 4) - 1
      do k = 0, ke_out - ks_out
        do j = 0, je_out - js_out
          do i = 0, ie_out - is_out
            if (data_out(is_out + i, js_out + j, ks_out + k, l + 1, 1) .lt. &
              data_in(is_in + i, js_in + j, ks_in + k, l + 1) ) then
                data_out(is_out + i, js_out + j, ks_out + k, l + 1, 1) = &
                  data_in(is_in +i, js_in + j, ks_in + k, l + 1)
            endif
          enddo
        enddo
      enddo
    enddo
  endif
end subroutine do_time_max_r4

!> Update the output buffer for reductions that involve summation (sum, avg, rms, pow).
!! Elements of the running field output buffer (data_out) are set with the following:
!!
!!    buffer(l) = buffer(l) + (weight * field(l)) ^ pow
!!
!! Where l are the indices passed in through the bounds_in/out
subroutine do_time_sum_update_r4(data_out, weight_sum, data_in, mask, is_masked, mask_variant, bounds_in, bounds_out, &
                               missing_value, diurnal_section, weight, pow)
  real(r4_kind),       intent(inout) :: data_out(:,:,:,:,:) !< output data
  real(r8_kind),             intent(inout) :: weight_sum(:,:,:,:) !< Sum of weights from the output buffer object
  real(r4_kind),       intent(in)    :: data_in(:,:,:,:)    !< data to update the buffer with
  logical,                   intent(in)    :: mask(:,:,:,:)       !< mask
  logical,                   intent(in)    :: is_masked           !< .True. if the field is using a mask
  logical,                   intent(in)    :: mask_variant        !< .True. if the mask changes over time
  type(fmsDiagIbounds_type), intent(in)    :: bounds_in           !< indices indicating the correct portion
                                                                  !! of the input buffer
  type(fmsDiagIbounds_type), intent(in)    :: bounds_out          !< indices indicating the correct portion
                                                                  !! of the output buffer
  real(r4_kind),       intent(in)    :: missing_value       !< Missing_value for data points that are masked
  integer, intent(in)                      :: diurnal_section !< the diurnal "section" if doing a diurnal reduction
                                                              !! indicates which index to add data on 5th axis
                                                              !! if not doing a diurnal reduction, this should always =1
  real(r8_kind),optional, intent(in)       :: weight          !< Weight applied to data_in before added to data_out
                                                              !! used for weighted averages, default 1.0
  integer ,optional, intent(in) :: pow                            !< Used for pow(er) reduction,
                                                                  !! calculates field_data^pow before adding to buffer

  real(r4_kind) :: weight_scale !< local copy of optional weight
  integer, parameter  :: kindl = r4_kind !< real kind size as set by macro
  integer :: diurnal !< diurnal index to indicate which daily section is updated
                     !! will be 1 unless using a diurnal reduction

  if(present(weight)) then
    weight_scale = real(weight, kind=kindl)
  else
    weight_scale = 1.0_kindl
  endif

  if(diurnal_section .lt. 0) then
    diurnal = 1
  else
    diurnal = diurnal_section
  endif

  if (is_masked) then
    if (mask_variant) then
      ! Mask changes over time so the weight is an array
      call sum_mask_variant(data_out, data_in, weight_sum, bounds_in, bounds_out, mask, diurnal, weight_scale, pow)
    else
      call sum_mask(data_out, data_in, weight_sum, bounds_in, bounds_out, mask, diurnal, &
        missing_value, weight_scale, pow)
    endif
  else
    call sum_no_mask(data_out, data_in, weight_sum, bounds_in, bounds_out, diurnal, weight_scale, pow)
  endif
end subroutine do_time_sum_update_r4

subroutine sum_mask_r4(data_out, data_in, weight_sum, bounds_in, bounds_out, mask, diurnal, missing_value, &
  weight_scale, pow)
  real(r4_kind),       intent(inout) :: data_out(:,:,:,:,:) !< output data
  real(r4_kind),       intent(in)    :: data_in(:,:,:,:)    !< data to update the buffer with
  real(r8_kind),             intent(inout) :: weight_sum(:,:,:,:) !< Sum of weights from the output buffer object
  type(fmsDiagIbounds_type), intent(in)    :: bounds_in           !< indices indicating the correct portion
                                                                  !! of the input buffer
  type(fmsDiagIbounds_type), intent(in)    :: bounds_out          !< indices indicating the correct portion
                                                                  !! of the output buffer
  logical,                   intent(in)    :: mask(:,:,:,:)       !< mask
  integer,                   intent(in)    :: diurnal             !< diurnal index to indicate which daily section is
                                                                  !! updated will be 1 unless using a diurnal reduction
  real(r4_kind),       intent(in)    :: missing_value       !< Missing_value for data points that are masked
  real(r4_kind),       intent(in)    :: weight_scale        !< weight scale to use
  integer ,optional,         intent(in)    :: pow                 !< Used for pow(er) reduction,
                                                                  !! calculates field_data^pow before adding to buffer

  integer :: is_in, ie_in, js_in, je_in, ks_in, ke_in       !< Starting and ending indices of each dimention for
                                                            !! the input buffer
  integer :: is_out, ie_out, js_out, je_out, ks_out, ke_out !< Starting and ending indices of each dimention for
                                                            !! the output buffer
  integer :: pow_loc !> local copy of optional pow value (set if using pow reduction)
  integer :: i, j, k, l !< For looping

  is_out = bounds_out%get_imin()
  ie_out = bounds_out%get_imax()
  js_out = bounds_out%get_jmin()
  je_out = bounds_out%get_jmax()
  ks_out = bounds_out%get_kmin()
  ke_out = bounds_out%get_kmax()

  is_in = bounds_in%get_imin()
  ie_in = bounds_in%get_imax()
  js_in = bounds_in%get_jmin()
  je_in = bounds_in%get_jmax()
  ks_in = bounds_in%get_kmin()
  ke_in = bounds_in%get_kmax()

  weight_sum = weight_sum + weight_scale
  if (present(pow)) then
    do k = 0, ke_out - ks_out
      do j = 0, je_out - js_out
        do i = 0, ie_out - is_out
          where (mask(is_in + i, js_in + j, ks_in + k, :))
            data_out(is_out + i, js_out + j, ks_out + k, :, diurnal) =           &
                   data_out(is_out + i, js_out + j, ks_out + k, :, diurnal)  &
                 + (data_in(is_in +i, js_in + j, ks_in + k, :) * weight_scale) ** pow
          elsewhere
            data_out(is_out + i, js_out + j, ks_out + k, :, diurnal) = missing_value
          endwhere
        enddo
      enddo
    enddo
  else
    do k = 0, ke_out - ks_out
      do j = 0, je_out - js_out
        do i = 0, ie_out - is_out
          where (mask(is_in + i, js_in + j, ks_in + k, :))
            data_out(is_out + i, js_out + j, ks_out + k, :, diurnal) =           &
                   data_out(is_out + i, js_out + j, ks_out + k, :, diurnal)  &
                 + (data_in(is_in +i, js_in + j, ks_in + k, :) * weight_scale)
          elsewhere
            data_out(is_out + i, js_out + j, ks_out + k, :, diurnal) = missing_value
          endwhere
        enddo
      enddo
    enddo
  endif
end subroutine sum_mask_r4

subroutine sum_mask_variant_r4(data_out, data_in, weight_sum, bounds_in, bounds_out, mask, diurnal, weight_scale, pow)
  real(r4_kind),       intent(inout) :: data_out(:,:,:,:,:) !< output data
  real(r4_kind),       intent(in)    :: data_in(:,:,:,:)    !< data to update the buffer with
  real(r8_kind),             intent(inout) :: weight_sum(:,:,:,:) !< Sum of weights from the output buffer object
  type(fmsDiagIbounds_type), intent(in)    :: bounds_in           !< indices indicating the correct portion
                                                                  !! of the input buffer
  type(fmsDiagIbounds_type), intent(in)    :: bounds_out          !< indices indicating the correct portion
                                                                  !! of the output buffer
  logical,                   intent(in)    :: mask(:,:,:,:)       !< mask
  integer,                   intent(in)    :: diurnal             !< diurnal index to indicate which daily section is
                                                                  !! updated will be 1 unless using a diurnal reduction
  real(r4_kind),       intent(in)    :: weight_scale        !< weight scale to use
  integer ,optional,         intent(in)    :: pow                 !< Used for pow(er) reduction,
                                                                  !! calculates field_data^pow before adding to buffer

  integer :: is_in, ie_in, js_in, je_in, ks_in, ke_in       !< Starting and ending indices of each dimention for
                                                            !! the input buffer
  integer :: is_out, ie_out, js_out, je_out, ks_out, ke_out !< Starting and ending indices of each dimention for
                                                            !! the output buffer
  integer :: pow_loc !> local copy of optional pow value (set if using pow reduction)
  integer :: i, j, k, l !< For looping

  is_out = bounds_out%get_imin()
  ie_out = bounds_out%get_imax()
  js_out = bounds_out%get_jmin()
  je_out = bounds_out%get_jmax()
  ks_out = bounds_out%get_kmin()
  ke_out = bounds_out%get_kmax()

  is_in = bounds_in%get_imin()
  ie_in = bounds_in%get_imax()
  js_in = bounds_in%get_jmin()
  je_in = bounds_in%get_jmax()
  ks_in = bounds_in%get_kmin()
  ke_in = bounds_in%get_kmax()

  if (present(pow)) then
    do k = 0, ke_out - ks_out
      do j = 0, je_out - js_out
        do i = 0, ie_out - is_out
          where (mask(is_in + i, js_in + j, ks_in + k, :))
            data_out(is_out + i, js_out + j, ks_out + k, :, diurnal) =           &
                   data_out(is_out + i, js_out + j, ks_out + k, :, diurnal)  &
                 + (data_in(is_in +i, js_in + j, ks_in + k, :) * weight_scale) ** pow

            !Increase the weight sum for the grid point that was not masked
            weight_sum(is_out + i, js_out + j, ks_out + k, :) = &
              weight_sum(is_out + i, js_out + j, ks_out + k, :) + weight_scale
          endwhere
        enddo
      enddo
    enddo
  else
    do k = 0, ke_out - ks_out
      do j = 0, je_out - js_out
        do i = 0, ie_out - is_out
          where (mask(is_in + i, js_in + j, ks_in + k, :))
            data_out(is_out + i, js_out + j, ks_out + k, :, diurnal) =           &
                   data_out(is_out + i, js_out + j, ks_out + k, :, diurnal)  &
                 + (data_in(is_in +i, js_in + j, ks_in + k, :) * weight_scale)

            !Increase the weight sum for the grid point that was not masked
            weight_sum(is_out + i, js_out + j, ks_out + k, :) = &
              weight_sum(is_out + i, js_out + j, ks_out + k, :) + weight_scale
          endwhere
        enddo
      enddo
    enddo
  endif
end subroutine sum_mask_variant_r4

subroutine sum_no_mask_r4(data_out, data_in, weight_sum, bounds_in, bounds_out, diurnal, weight_scale, pow)
  real(r4_kind),       intent(inout) :: data_out(:,:,:,:,:) !< output data
  real(r4_kind),       intent(in)    :: data_in(:,:,:,:)    !< data to update the buffer with
  real(r8_kind),             intent(inout) :: weight_sum(:,:,:,:) !< Sum of weights from the output buffer object
  type(fmsDiagIbounds_type), intent(in)    :: bounds_in           !< indices indicating the correct portion
                                                                  !! of the input buffer
  type(fmsDiagIbounds_type), intent(in)    :: bounds_out          !< indices indicating the correct portion
                                                                  !! of the output buffer
  integer,                   intent(in)    :: diurnal             !< diurnal index to indicate which daily section is
                                                                  !! updated will be 1 unless using a diurnal reduction
  real(r4_kind),       intent(in)    :: weight_scale        !< weight scale to use
  integer ,optional,         intent(in)    :: pow                 !< Used for pow(er) reduction,
                                                                  !! calculates field_data^pow before adding to buffer

  integer :: is_in, ie_in, js_in, je_in, ks_in, ke_in       !< Starting and ending indices of each dimention for
                                                            !! the input buffer
  integer :: is_out, ie_out, js_out, je_out, ks_out, ke_out !< Starting and ending indices of each dimention for
                                                            !! the output buffer
  integer :: i, j, k, l !< For looping

  is_out = bounds_out%get_imin()
  ie_out = bounds_out%get_imax()
  js_out = bounds_out%get_jmin()
  je_out = bounds_out%get_jmax()
  ks_out = bounds_out%get_kmin()
  ke_out = bounds_out%get_kmax()

  is_in = bounds_in%get_imin()
  ie_in = bounds_in%get_imax()
  js_in = bounds_in%get_jmin()
  je_in = bounds_in%get_jmax()
  ks_in = bounds_in%get_kmin()
  ke_in = bounds_in%get_kmax()

  weight_sum = weight_sum + weight_scale

  if (present(pow)) then
    do k = 0, ke_out - ks_out
      do j = 0, je_out - js_out
        do i = 0, ie_out - is_out
          data_out(is_out + i, js_out + j, ks_out + k,  :, diurnal) =  &
            data_out(is_out + i, js_out + j, ks_out + k, :, diurnal)  &
            + (data_in(is_in +i, js_in + j, ks_in + k, :) * weight_scale) ** pow
        enddo
      enddo
    enddo
  else
    do k = 0, ke_out - ks_out
      do j = 0, je_out - js_out
        do i = 0, ie_out - is_out
          data_out(is_out + i, js_out + j, ks_out + k,  :, diurnal) =  &
            data_out(is_out + i, js_out + j, ks_out + k, :, diurnal)  &
            + (data_in(is_in +i, js_in + j, ks_in + k, :) * weight_scale)
        enddo
      enddo
    enddo
  endif
end subroutine sum_no_mask_r4

!> To be called with diag_send_complete, finishes reductions
!! Just divides the buffer by the counter array(which is just the sum of the weights used in the buffer's reduction)
!! TODO: change has_mask to an actual logical mask so we don't have to check for missing values
subroutine sum_update_done_r4(out_buffer_data, weight_sum, reduction_method, missing_val, has_mask, mask_variant, &
                            n_diurnal_samples)
  real(r4_kind), intent(inout) :: out_buffer_data(:,:,:,:,:) !< data buffer previously updated with
                                                                   !! do_time_sum_update
  real(r8_kind), intent(in)          :: weight_sum(:,:,:,:) !< sum of weights for averaging,
                                                            !! provided via argument to send data
  integer, intent(in)                :: reduction_method !< which reduction method to use
                                                         !! should always be one of time_avg, time_diurnal, or time_rms
  real(r4_kind), intent(in)    :: missing_val !< missing value for masked elements
  logical, intent(in)                :: has_mask !< indicates if mask is used so missing values can be skipped
  logical, intent(in)                :: mask_variant !< Indicates if the mask changes over time
  integer, optional, intent(in)      :: n_diurnal_samples !< number of diurnal samples as set in reduction method
  integer, allocatable :: wsum(:,:,:,:) !< local cp of weight_sum, only changed if using diurnal
  !! TODO replace conditional in the `where` with passed in and ajusted mask from the original call
  !logical, optional, intent(in)      :: mask(:,:,:,:) !< logical mask from accept data call, if using one.
  !logical                            :: has_mask !< whether or not mask is present

  integer :: i, j, k, l !< For do loops

  allocate(wsum(size(weight_sum,1), size(weight_sum,3), size(weight_sum,3), size(weight_sum,4)))
  ! need to divide weight sum by amount of samples to get the actual
  ! number of times that the diurnal section was incremented
  ! legacy diag manager stored these weights explicitly, this doesn't so assumes uniformity in when data is sent
  if(reduction_method .eq. time_diurnal) then
    if(.not. present(n_diurnal_samples)) call mpp_error(FATAL, &
      "SUM_UPDATE_DONE_ :: reduction method is diurnal but no sample size was given")
    wsum = weight_sum / n_diurnal_samples
  else
    wsum = weight_sum
  endif

  if ( has_mask ) then
    if (.not. mask_variant) then
      ! The mask does not change over time so wsum is just an integer and it is the same value for all fields
      where(out_buffer_data(:,:,:,:,:) .ne. missing_val)
        out_buffer_data(:,:,:,:,:) = out_buffer_data(:,:,:,:,:) &
                                   / wsum(1,1,1,1)
      endwhere
    else
      ! The mask changes over time
      do l = 1, size(out_buffer_data, 4)
        do k = 1, size(out_buffer_data, 3)
          do j = 1, size(out_buffer_data, 2)
            do i = 1, size(out_buffer_data, 1)
              if (wsum(i, j, k, l) .gt. 0) then
                out_buffer_data(i,j,k,l,:) = out_buffer_data(i,j,k,l,:)/ wsum(i,j,k,l)
              else
                ! Data was never received
                out_buffer_data(i,j,k,l,:) = missing_val
              endif
            enddo
          enddo
        enddo
      enddo
    endif
  else
    ! There is no mask!
    out_buffer_data(:,:,:,:,:) = out_buffer_data(:,:,:,:,:) &
                               / wsum(1,1,1,1)
  endif

  if(reduction_method .eq. time_rms .and. has_mask) then
    where(out_buffer_data(:,:,:,:,1) .ne. missing_val)
      out_buffer_data(:,:,:,:,1) = SQRT(out_buffer_data(:,:,:,:,1))
    endwhere
  else if(reduction_method .eq. time_rms) then
    out_buffer_data(:,:,:,:,1) = SQRT(out_buffer_data(:,:,:,:,1))
  endif

end subroutine

# 47 "diag_manager/include/fms_diag_reduction_methods_r4.fh" 2
# 181 "diag_manager/fms_diag_reduction_methods.F90" 2

# 1 "diag_manager/include/fms_diag_reduction_methods_r8.fh" 1
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





























# 1 "diag_manager/include/fms_diag_reduction_methods.inc" 1
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

! for any debug prints




!> @brief Do the time_none reduction method (i.e copy the correct portion of the input data)
subroutine do_time_none_r8 (data_out, data_in, mask, is_masked, bounds_in, bounds_out, missing_value)
  real(r8_kind),       intent(inout) :: data_out(:,:,:,:,:) !< output data
  real(r8_kind),       intent(in)    :: data_in(:,:,:,:)    !< data to update the buffer with
  logical,                   intent(in)    :: mask(:,:,:,:)       !< mask
  logical,                   intent(in)    :: is_masked           !< .True. if the field is using a mask
  type(fmsDiagIbounds_type), intent(in)    :: bounds_in           !< indices indicating the correct portion
                                                                  !! of the input buffer
  type(fmsDiagIbounds_type), intent(in)    :: bounds_out          !< indices indicating the correct portion
                                                                  !! of the output buffer
  real(r8_kind),       intent(in)    :: missing_value       !< Missing_value for data points that are masked

  integer :: is_in, ie_in, js_in, je_in, ks_in, ke_in       !< Starting and ending indices of each dimention for
                                                            !! the input buffer
  integer :: is_out, ie_out, js_out, je_out, ks_out, ke_out !< Starting and ending indices of each dimention for
                                                            !! the output buffer

  is_out = bounds_out%get_imin()
  ie_out = bounds_out%get_imax()
  js_out = bounds_out%get_jmin()
  je_out = bounds_out%get_jmax()
  ks_out = bounds_out%get_kmin()
  ke_out = bounds_out%get_kmax()

  is_in = bounds_in%get_imin()
  ie_in = bounds_in%get_imax()
  js_in = bounds_in%get_jmin()
  je_in = bounds_in%get_jmax()
  ks_in = bounds_in%get_kmin()
  ke_in = bounds_in%get_kmax()

  if (is_masked) then
    where (mask(is_in:ie_in, js_in:je_in, ks_in:ke_in, :))
      data_out(is_out:ie_out, js_out:je_out, ks_out:ke_out, :, 1) = &
      data_in(is_in:ie_in, js_in:je_in, ks_in:ke_in, :)
    elsewhere
      data_out(is_out:ie_out, js_out:je_out, ks_out:ke_out, :, 1) = missing_value
    end where
  else
    data_out(is_out:ie_out, js_out:je_out, ks_out:ke_out, :, 1) = &
      data_in(is_in:ie_in, js_in:je_in, ks_in:ke_in, :)
  endif

end subroutine do_time_none_r8

!> @brief Do the time_min reduction method (i.e maintain the minimum value of the averaging time)
subroutine do_time_min_r8 (data_out, data_in, mask, is_masked, bounds_in, bounds_out, missing_value)
  real(r8_kind),       intent(inout) :: data_out(:,:,:,:,:) !< output data
  real(r8_kind),       intent(in)    :: data_in(:,:,:,:)    !< data to update the buffer with
  logical,                   intent(in)    :: mask(:,:,:,:)       !< mask
  logical,                   intent(in)    :: is_masked           !< .True. if the field is using a mask
  type(fmsDiagIbounds_type), intent(in)    :: bounds_in           !< indices indicating the correct portion
                                                                  !! of the input buffer
  type(fmsDiagIbounds_type), intent(in)    :: bounds_out          !< indices indicating the correct portion
                                                                  !! of the output buffer
  real(r8_kind),       intent(in)    :: missing_value       !< Missing_value for data points that are masked

  integer :: is_in, ie_in, js_in, je_in, ks_in, ke_in       !< Starting and ending indices of each dimention for
                                                            !! the input buffer
  integer :: is_out, ie_out, js_out, je_out, ks_out, ke_out !< Starting and ending indices of each dimention for
                                                            !! the output buffer

  integer :: i, j, k, l !< For looping

  is_out = bounds_out%get_imin()
  ie_out = bounds_out%get_imax()
  js_out = bounds_out%get_jmin()
  je_out = bounds_out%get_jmax()
  ks_out = bounds_out%get_kmin()
  ke_out = bounds_out%get_kmax()

  is_in = bounds_in%get_imin()
  ie_in = bounds_in%get_imax()
  js_in = bounds_in%get_jmin()
  je_in = bounds_in%get_jmax()
  ks_in = bounds_in%get_kmin()
  ke_in = bounds_in%get_kmax()

  !> Separated this loops for performance. If is_masked = .false. (i.e "mask" and "rmask" were never passed in)
  !! then mask will always be .True. so the if (mask) is redudant.
  if (is_masked) then
    do l = 0, size(data_out, 4) - 1
      do k = 0, ke_out - ks_out
        do j = 0, je_out - js_out
          do i = 0, ie_out - is_out
            if (mask(is_in + i, js_in + j, ks_in + k, l + 1)) then
              if (data_out(is_out + i, js_out + j, ks_out + k, l + 1, 1) .gt. &
                data_in(is_in + i, js_in + j, ks_in + k, l + 1) ) then
                  data_out(is_out + i, js_out + j, ks_out + k, l + 1, 1) = &
                    data_in(is_in +i, js_in + j, ks_in + k, l + 1)
              endif
            else
              data_out(is_out + i, js_out + j, ks_out + k, l + 1, 1) = missing_value
            endif
          enddo
        enddo
      enddo
    enddo
  else
    do l = 0, size(data_out, 4) - 1
      do k = 0, ke_out - ks_out
        do j = 0, je_out - js_out
          do i = 0, ie_out - is_out
            if (data_out(is_out + i, js_out + j, ks_out + k, l + 1, 1) .gt. &
              data_in(is_in + i, js_in + j, ks_in + k, l + 1) ) then
                data_out(is_out + i, js_out + j, ks_out + k, l + 1, 1) = &
                  data_in(is_in +i, js_in + j, ks_in + k, l + 1)
            endif
          enddo
        enddo
      enddo
    enddo
  endif

end subroutine do_time_min_r8

!> @brief Do the time_max reduction method (i.e maintain the maximum value of the averaging time)
subroutine do_time_max_r8 (data_out, data_in, mask, is_masked, bounds_in, bounds_out, missing_value)
  real(r8_kind),       intent(inout) :: data_out(:,:,:,:,:) !< output data
  real(r8_kind),       intent(in)    :: data_in(:,:,:,:)    !< data to update the buffer with
  logical,                   intent(in)    :: mask(:,:,:,:)       !< mask
  logical,                   intent(in)    :: is_masked           !< .True. if the field is using a mask
  type(fmsDiagIbounds_type), intent(in)    :: bounds_in           !< indices indicating the correct portion
                                                                  !! of the input buffer
  type(fmsDiagIbounds_type), intent(in)    :: bounds_out          !< indices indicating the correct portion
                                                                  !! of the output buffer
  real(r8_kind),       intent(in)    :: missing_value       !< Missing_value for data points that are masked

  integer :: is_in, ie_in, js_in, je_in, ks_in, ke_in       !< Starting and ending indices of each dimention for
                                                            !! the input buffer
  integer :: is_out, ie_out, js_out, je_out, ks_out, ke_out !< Starting and ending indices of each dimention for
                                                            !! the output buffer

  integer :: i, j, k, l !< For looping

  is_out = bounds_out%get_imin()
  ie_out = bounds_out%get_imax()
  js_out = bounds_out%get_jmin()
  je_out = bounds_out%get_jmax()
  ks_out = bounds_out%get_kmin()
  ke_out = bounds_out%get_kmax()

  is_in = bounds_in%get_imin()
  ie_in = bounds_in%get_imax()
  js_in = bounds_in%get_jmin()
  je_in = bounds_in%get_jmax()
  ks_in = bounds_in%get_kmin()
  ke_in = bounds_in%get_kmax()

  !> Separated this loops for performance. If is_masked = .false. (i.e "mask" and "rmask" were never passed in)
  !! then mask will always be .True. so the if (mask) is redudant.
  if (is_masked) then
    do l = 0, size(data_out, 4) - 1
      do k = 0, ke_out - ks_out
        do j = 0, je_out - js_out
          do i = 0, ie_out - is_out
            if (mask(is_in + i, js_in + j, ks_in + k, l + 1)) then
              if (data_out(is_out + i, js_out + j, ks_out + k, l + 1, 1) .lt. &
                data_in(is_in + i, js_in + j, ks_in + k, l + 1) ) then
                  data_out(is_out + i, js_out + j, ks_out + k, l + 1, 1) = &
                    data_in(is_in +i, js_in + j, ks_in + k, l + 1)
              endif
            else
              data_out(is_out + i, js_out + j, ks_out + k, l + 1, 1) = missing_value
            endif
          enddo
        enddo
      enddo
    enddo
  else
    do l = 0, size(data_out, 4) - 1
      do k = 0, ke_out - ks_out
        do j = 0, je_out - js_out
          do i = 0, ie_out - is_out
            if (data_out(is_out + i, js_out + j, ks_out + k, l + 1, 1) .lt. &
              data_in(is_in + i, js_in + j, ks_in + k, l + 1) ) then
                data_out(is_out + i, js_out + j, ks_out + k, l + 1, 1) = &
                  data_in(is_in +i, js_in + j, ks_in + k, l + 1)
            endif
          enddo
        enddo
      enddo
    enddo
  endif
end subroutine do_time_max_r8

!> Update the output buffer for reductions that involve summation (sum, avg, rms, pow).
!! Elements of the running field output buffer (data_out) are set with the following:
!!
!!    buffer(l) = buffer(l) + (weight * field(l)) ^ pow
!!
!! Where l are the indices passed in through the bounds_in/out
subroutine do_time_sum_update_r8(data_out, weight_sum, data_in, mask, is_masked, mask_variant, bounds_in, bounds_out, &
                               missing_value, diurnal_section, weight, pow)
  real(r8_kind),       intent(inout) :: data_out(:,:,:,:,:) !< output data
  real(r8_kind),             intent(inout) :: weight_sum(:,:,:,:) !< Sum of weights from the output buffer object
  real(r8_kind),       intent(in)    :: data_in(:,:,:,:)    !< data to update the buffer with
  logical,                   intent(in)    :: mask(:,:,:,:)       !< mask
  logical,                   intent(in)    :: is_masked           !< .True. if the field is using a mask
  logical,                   intent(in)    :: mask_variant        !< .True. if the mask changes over time
  type(fmsDiagIbounds_type), intent(in)    :: bounds_in           !< indices indicating the correct portion
                                                                  !! of the input buffer
  type(fmsDiagIbounds_type), intent(in)    :: bounds_out          !< indices indicating the correct portion
                                                                  !! of the output buffer
  real(r8_kind),       intent(in)    :: missing_value       !< Missing_value for data points that are masked
  integer, intent(in)                      :: diurnal_section !< the diurnal "section" if doing a diurnal reduction
                                                              !! indicates which index to add data on 5th axis
                                                              !! if not doing a diurnal reduction, this should always =1
  real(r8_kind),optional, intent(in)       :: weight          !< Weight applied to data_in before added to data_out
                                                              !! used for weighted averages, default 1.0
  integer ,optional, intent(in) :: pow                            !< Used for pow(er) reduction,
                                                                  !! calculates field_data^pow before adding to buffer

  real(r8_kind) :: weight_scale !< local copy of optional weight
  integer, parameter  :: kindl = r8_kind !< real kind size as set by macro
  integer :: diurnal !< diurnal index to indicate which daily section is updated
                     !! will be 1 unless using a diurnal reduction

  if(present(weight)) then
    weight_scale = real(weight, kind=kindl)
  else
    weight_scale = 1.0_kindl
  endif

  if(diurnal_section .lt. 0) then
    diurnal = 1
  else
    diurnal = diurnal_section
  endif

  if (is_masked) then
    if (mask_variant) then
      ! Mask changes over time so the weight is an array
      call sum_mask_variant(data_out, data_in, weight_sum, bounds_in, bounds_out, mask, diurnal, weight_scale, pow)
    else
      call sum_mask(data_out, data_in, weight_sum, bounds_in, bounds_out, mask, diurnal, &
        missing_value, weight_scale, pow)
    endif
  else
    call sum_no_mask(data_out, data_in, weight_sum, bounds_in, bounds_out, diurnal, weight_scale, pow)
  endif
end subroutine do_time_sum_update_r8

subroutine sum_mask_r8(data_out, data_in, weight_sum, bounds_in, bounds_out, mask, diurnal, missing_value, &
  weight_scale, pow)
  real(r8_kind),       intent(inout) :: data_out(:,:,:,:,:) !< output data
  real(r8_kind),       intent(in)    :: data_in(:,:,:,:)    !< data to update the buffer with
  real(r8_kind),             intent(inout) :: weight_sum(:,:,:,:) !< Sum of weights from the output buffer object
  type(fmsDiagIbounds_type), intent(in)    :: bounds_in           !< indices indicating the correct portion
                                                                  !! of the input buffer
  type(fmsDiagIbounds_type), intent(in)    :: bounds_out          !< indices indicating the correct portion
                                                                  !! of the output buffer
  logical,                   intent(in)    :: mask(:,:,:,:)       !< mask
  integer,                   intent(in)    :: diurnal             !< diurnal index to indicate which daily section is
                                                                  !! updated will be 1 unless using a diurnal reduction
  real(r8_kind),       intent(in)    :: missing_value       !< Missing_value for data points that are masked
  real(r8_kind),       intent(in)    :: weight_scale        !< weight scale to use
  integer ,optional,         intent(in)    :: pow                 !< Used for pow(er) reduction,
                                                                  !! calculates field_data^pow before adding to buffer

  integer :: is_in, ie_in, js_in, je_in, ks_in, ke_in       !< Starting and ending indices of each dimention for
                                                            !! the input buffer
  integer :: is_out, ie_out, js_out, je_out, ks_out, ke_out !< Starting and ending indices of each dimention for
                                                            !! the output buffer
  integer :: pow_loc !> local copy of optional pow value (set if using pow reduction)
  integer :: i, j, k, l !< For looping

  is_out = bounds_out%get_imin()
  ie_out = bounds_out%get_imax()
  js_out = bounds_out%get_jmin()
  je_out = bounds_out%get_jmax()
  ks_out = bounds_out%get_kmin()
  ke_out = bounds_out%get_kmax()

  is_in = bounds_in%get_imin()
  ie_in = bounds_in%get_imax()
  js_in = bounds_in%get_jmin()
  je_in = bounds_in%get_jmax()
  ks_in = bounds_in%get_kmin()
  ke_in = bounds_in%get_kmax()

  weight_sum = weight_sum + weight_scale
  if (present(pow)) then
    do k = 0, ke_out - ks_out
      do j = 0, je_out - js_out
        do i = 0, ie_out - is_out
          where (mask(is_in + i, js_in + j, ks_in + k, :))
            data_out(is_out + i, js_out + j, ks_out + k, :, diurnal) =           &
                   data_out(is_out + i, js_out + j, ks_out + k, :, diurnal)  &
                 + (data_in(is_in +i, js_in + j, ks_in + k, :) * weight_scale) ** pow
          elsewhere
            data_out(is_out + i, js_out + j, ks_out + k, :, diurnal) = missing_value
          endwhere
        enddo
      enddo
    enddo
  else
    do k = 0, ke_out - ks_out
      do j = 0, je_out - js_out
        do i = 0, ie_out - is_out
          where (mask(is_in + i, js_in + j, ks_in + k, :))
            data_out(is_out + i, js_out + j, ks_out + k, :, diurnal) =           &
                   data_out(is_out + i, js_out + j, ks_out + k, :, diurnal)  &
                 + (data_in(is_in +i, js_in + j, ks_in + k, :) * weight_scale)
          elsewhere
            data_out(is_out + i, js_out + j, ks_out + k, :, diurnal) = missing_value
          endwhere
        enddo
      enddo
    enddo
  endif
end subroutine sum_mask_r8

subroutine sum_mask_variant_r8(data_out, data_in, weight_sum, bounds_in, bounds_out, mask, diurnal, weight_scale, pow)
  real(r8_kind),       intent(inout) :: data_out(:,:,:,:,:) !< output data
  real(r8_kind),       intent(in)    :: data_in(:,:,:,:)    !< data to update the buffer with
  real(r8_kind),             intent(inout) :: weight_sum(:,:,:,:) !< Sum of weights from the output buffer object
  type(fmsDiagIbounds_type), intent(in)    :: bounds_in           !< indices indicating the correct portion
                                                                  !! of the input buffer
  type(fmsDiagIbounds_type), intent(in)    :: bounds_out          !< indices indicating the correct portion
                                                                  !! of the output buffer
  logical,                   intent(in)    :: mask(:,:,:,:)       !< mask
  integer,                   intent(in)    :: diurnal             !< diurnal index to indicate which daily section is
                                                                  !! updated will be 1 unless using a diurnal reduction
  real(r8_kind),       intent(in)    :: weight_scale        !< weight scale to use
  integer ,optional,         intent(in)    :: pow                 !< Used for pow(er) reduction,
                                                                  !! calculates field_data^pow before adding to buffer

  integer :: is_in, ie_in, js_in, je_in, ks_in, ke_in       !< Starting and ending indices of each dimention for
                                                            !! the input buffer
  integer :: is_out, ie_out, js_out, je_out, ks_out, ke_out !< Starting and ending indices of each dimention for
                                                            !! the output buffer
  integer :: pow_loc !> local copy of optional pow value (set if using pow reduction)
  integer :: i, j, k, l !< For looping

  is_out = bounds_out%get_imin()
  ie_out = bounds_out%get_imax()
  js_out = bounds_out%get_jmin()
  je_out = bounds_out%get_jmax()
  ks_out = bounds_out%get_kmin()
  ke_out = bounds_out%get_kmax()

  is_in = bounds_in%get_imin()
  ie_in = bounds_in%get_imax()
  js_in = bounds_in%get_jmin()
  je_in = bounds_in%get_jmax()
  ks_in = bounds_in%get_kmin()
  ke_in = bounds_in%get_kmax()

  if (present(pow)) then
    do k = 0, ke_out - ks_out
      do j = 0, je_out - js_out
        do i = 0, ie_out - is_out
          where (mask(is_in + i, js_in + j, ks_in + k, :))
            data_out(is_out + i, js_out + j, ks_out + k, :, diurnal) =           &
                   data_out(is_out + i, js_out + j, ks_out + k, :, diurnal)  &
                 + (data_in(is_in +i, js_in + j, ks_in + k, :) * weight_scale) ** pow

            !Increase the weight sum for the grid point that was not masked
            weight_sum(is_out + i, js_out + j, ks_out + k, :) = &
              weight_sum(is_out + i, js_out + j, ks_out + k, :) + weight_scale
          endwhere
        enddo
      enddo
    enddo
  else
    do k = 0, ke_out - ks_out
      do j = 0, je_out - js_out
        do i = 0, ie_out - is_out
          where (mask(is_in + i, js_in + j, ks_in + k, :))
            data_out(is_out + i, js_out + j, ks_out + k, :, diurnal) =           &
                   data_out(is_out + i, js_out + j, ks_out + k, :, diurnal)  &
                 + (data_in(is_in +i, js_in + j, ks_in + k, :) * weight_scale)

            !Increase the weight sum for the grid point that was not masked
            weight_sum(is_out + i, js_out + j, ks_out + k, :) = &
              weight_sum(is_out + i, js_out + j, ks_out + k, :) + weight_scale
          endwhere
        enddo
      enddo
    enddo
  endif
end subroutine sum_mask_variant_r8

subroutine sum_no_mask_r8(data_out, data_in, weight_sum, bounds_in, bounds_out, diurnal, weight_scale, pow)
  real(r8_kind),       intent(inout) :: data_out(:,:,:,:,:) !< output data
  real(r8_kind),       intent(in)    :: data_in(:,:,:,:)    !< data to update the buffer with
  real(r8_kind),             intent(inout) :: weight_sum(:,:,:,:) !< Sum of weights from the output buffer object
  type(fmsDiagIbounds_type), intent(in)    :: bounds_in           !< indices indicating the correct portion
                                                                  !! of the input buffer
  type(fmsDiagIbounds_type), intent(in)    :: bounds_out          !< indices indicating the correct portion
                                                                  !! of the output buffer
  integer,                   intent(in)    :: diurnal             !< diurnal index to indicate which daily section is
                                                                  !! updated will be 1 unless using a diurnal reduction
  real(r8_kind),       intent(in)    :: weight_scale        !< weight scale to use
  integer ,optional,         intent(in)    :: pow                 !< Used for pow(er) reduction,
                                                                  !! calculates field_data^pow before adding to buffer

  integer :: is_in, ie_in, js_in, je_in, ks_in, ke_in       !< Starting and ending indices of each dimention for
                                                            !! the input buffer
  integer :: is_out, ie_out, js_out, je_out, ks_out, ke_out !< Starting and ending indices of each dimention for
                                                            !! the output buffer
  integer :: i, j, k, l !< For looping

  is_out = bounds_out%get_imin()
  ie_out = bounds_out%get_imax()
  js_out = bounds_out%get_jmin()
  je_out = bounds_out%get_jmax()
  ks_out = bounds_out%get_kmin()
  ke_out = bounds_out%get_kmax()

  is_in = bounds_in%get_imin()
  ie_in = bounds_in%get_imax()
  js_in = bounds_in%get_jmin()
  je_in = bounds_in%get_jmax()
  ks_in = bounds_in%get_kmin()
  ke_in = bounds_in%get_kmax()

  weight_sum = weight_sum + weight_scale

  if (present(pow)) then
    do k = 0, ke_out - ks_out
      do j = 0, je_out - js_out
        do i = 0, ie_out - is_out
          data_out(is_out + i, js_out + j, ks_out + k,  :, diurnal) =  &
            data_out(is_out + i, js_out + j, ks_out + k, :, diurnal)  &
            + (data_in(is_in +i, js_in + j, ks_in + k, :) * weight_scale) ** pow
        enddo
      enddo
    enddo
  else
    do k = 0, ke_out - ks_out
      do j = 0, je_out - js_out
        do i = 0, ie_out - is_out
          data_out(is_out + i, js_out + j, ks_out + k,  :, diurnal) =  &
            data_out(is_out + i, js_out + j, ks_out + k, :, diurnal)  &
            + (data_in(is_in +i, js_in + j, ks_in + k, :) * weight_scale)
        enddo
      enddo
    enddo
  endif
end subroutine sum_no_mask_r8

!> To be called with diag_send_complete, finishes reductions
!! Just divides the buffer by the counter array(which is just the sum of the weights used in the buffer's reduction)
!! TODO: change has_mask to an actual logical mask so we don't have to check for missing values
subroutine sum_update_done_r8(out_buffer_data, weight_sum, reduction_method, missing_val, has_mask, mask_variant, &
                            n_diurnal_samples)
  real(r8_kind), intent(inout) :: out_buffer_data(:,:,:,:,:) !< data buffer previously updated with
                                                                   !! do_time_sum_update
  real(r8_kind), intent(in)          :: weight_sum(:,:,:,:) !< sum of weights for averaging,
                                                            !! provided via argument to send data
  integer, intent(in)                :: reduction_method !< which reduction method to use
                                                         !! should always be one of time_avg, time_diurnal, or time_rms
  real(r8_kind), intent(in)    :: missing_val !< missing value for masked elements
  logical, intent(in)                :: has_mask !< indicates if mask is used so missing values can be skipped
  logical, intent(in)                :: mask_variant !< Indicates if the mask changes over time
  integer, optional, intent(in)      :: n_diurnal_samples !< number of diurnal samples as set in reduction method
  integer, allocatable :: wsum(:,:,:,:) !< local cp of weight_sum, only changed if using diurnal
  !! TODO replace conditional in the `where` with passed in and ajusted mask from the original call
  !logical, optional, intent(in)      :: mask(:,:,:,:) !< logical mask from accept data call, if using one.
  !logical                            :: has_mask !< whether or not mask is present

  integer :: i, j, k, l !< For do loops

  allocate(wsum(size(weight_sum,1), size(weight_sum,3), size(weight_sum,3), size(weight_sum,4)))
  ! need to divide weight sum by amount of samples to get the actual
  ! number of times that the diurnal section was incremented
  ! legacy diag manager stored these weights explicitly, this doesn't so assumes uniformity in when data is sent
  if(reduction_method .eq. time_diurnal) then
    if(.not. present(n_diurnal_samples)) call mpp_error(FATAL, &
      "SUM_UPDATE_DONE_ :: reduction method is diurnal but no sample size was given")
    wsum = weight_sum / n_diurnal_samples
  else
    wsum = weight_sum
  endif

  if ( has_mask ) then
    if (.not. mask_variant) then
      ! The mask does not change over time so wsum is just an integer and it is the same value for all fields
      where(out_buffer_data(:,:,:,:,:) .ne. missing_val)
        out_buffer_data(:,:,:,:,:) = out_buffer_data(:,:,:,:,:) &
                                   / wsum(1,1,1,1)
      endwhere
    else
      ! The mask changes over time
      do l = 1, size(out_buffer_data, 4)
        do k = 1, size(out_buffer_data, 3)
          do j = 1, size(out_buffer_data, 2)
            do i = 1, size(out_buffer_data, 1)
              if (wsum(i, j, k, l) .gt. 0) then
                out_buffer_data(i,j,k,l,:) = out_buffer_data(i,j,k,l,:)/ wsum(i,j,k,l)
              else
                ! Data was never received
                out_buffer_data(i,j,k,l,:) = missing_val
              endif
            enddo
          enddo
        enddo
      enddo
    endif
  else
    ! There is no mask!
    out_buffer_data(:,:,:,:,:) = out_buffer_data(:,:,:,:,:) &
                               / wsum(1,1,1,1)
  endif

  if(reduction_method .eq. time_rms .and. has_mask) then
    where(out_buffer_data(:,:,:,:,1) .ne. missing_val)
      out_buffer_data(:,:,:,:,1) = SQRT(out_buffer_data(:,:,:,:,1))
    endwhere
  else if(reduction_method .eq. time_rms) then
    out_buffer_data(:,:,:,:,1) = SQRT(out_buffer_data(:,:,:,:,1))
  endif

end subroutine

# 47 "diag_manager/include/fms_diag_reduction_methods_r8.fh" 2
# 182 "diag_manager/fms_diag_reduction_methods.F90" 2

end module fms_diag_reduction_methods_mod
!> @}
! close documentation grouping
