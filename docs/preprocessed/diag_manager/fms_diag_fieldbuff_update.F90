# 1 "diag_manager/fms_diag_fieldbuff_update.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "diag_manager/fms_diag_fieldbuff_update.F90"
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

!> @defgroup fms_diag_fieldbuff_update_mod fms_diag_fieldbuff_update_mod
!! @ingroup diag_manager
!! @{
!! @brief fms_diag_fieldbuff_update_mod Contains routines for updating the
!! buffer (array) of field data statistics (e.g. average, rms) with new field data.
!!
!!
!! <TT>fms_diag_fieldbuff_update_mod</TT> contains routines for updating the buffer
!!(array) of field data statistics (e.g. average, rms) with new field data. These
!! routines are called by the send_data routines in the diag_manager.
MODULE fms_diag_fieldbuff_update_mod
   USE platform_mod
   USE mpp_mod, ONLY: mpp_pe, mpp_root_pe
   USE time_manager_mod, ONLY: time_type
   USE fms_mod, ONLY: error_mesg, FATAL, WARNING, NOTE, stdout, stdlog, write_version_number,fms_error_handler
   USE diag_data_mod, ONLY:  debug_diag_manager
   USE fms_diag_outfield_mod, ONLY: fmsDiagOutfieldIndex_type, fmsDiagOutfield_type
   USE diag_util_mod, ONLY: fms_diag_check_out_of_bounds
   USE fms_diag_time_reduction_mod, ONLY: fmsDiagTimeReduction_type
   USE fms_diag_elem_weight_procs_mod, ONLY: addwf
   USE fms_diag_bbox_mod, ONLY: fmsDiagIbounds_type

   implicit none

   !!TODO: (MDM) Remove commented integer versions.

   !> @brief Interface fieldbuff_update updates elements of field output buffer based on input field
   !! data and mathematical operations on the field data.
   interface fieldbuff_update
      !< r4 version of the interface
      module procedure fieldbuff_update_r4
      !< r8 version of the interface
      module procedure fieldbuff_update_r8
      !< r4 version of the interface, where the field is 3D
      module procedure fieldbuff_update_3d_r4
      !< r8 version of the interface
      module procedure fieldbuff_update_3d_r8
      !< i4 version of the interface, , where the field is 3D
      !module procedure fieldbuff_update_i4
      !< i8 version of the interface
     ! module procedure fieldbuff_update_i8
   end interface

   !> @brief Interface fieldbuff_copy_missvals updates elements of the field output buffer with
   !! the missvalue input argument.
   interface fieldbuff_copy_missvals
      !< r4 version of the interface
      module procedure fieldbuff_copy_missvals_r4
      !< r8 version of the interface
      module procedure fieldbuff_copy_missvals_r8
      !< r4 version of the interface, , where the field is 3D
      module procedure fieldbuff_copy_missvals_3d_r4
      !< r8 version of the interface, , where the field is 3D
      module procedure fieldbuff_copy_missvals_3d_r8
      !< i4 version of the interface
      !module procedure fieldbuff_copy_missvals_i4
      !< i8 version of the interface
      !module procedure fieldbuff_copy_missvals_i8
   end interface

   !> @brief Interface fieldbuff_copy_fieldvals updates elements of the field output buffer with
   !! copies of corresponding element values in the input field data.
   interface fieldbuff_copy_fieldvals
      !< r4 version of the interface
      module procedure fieldbuff_copy_fieldvals_r4
      !< r8 version of the interface
      module procedure fieldbuff_copy_fieldvals_r8
      !< r4 version of the interface, , where the field is 3D
      module procedure fieldbuff_copy_fieldvals_3d_r4
      !< r8 version of the interface, , where the field is 3D
      module procedure fieldbuff_copy_fieldvals_3d_r8
      !< i4 version of the interface
      !module procedure fieldbuff_copy_fieldvals_i4
      !< i8 version of the interface
      !module procedure fieldbuff_copy_fieldvals_i8
  end interface
contains


# 1 "diag_manager/include/fms_diag_fieldbuff_update.inc" 1
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

# 34 "diag_manager/include/fms_diag_fieldbuff_update.inc"
# 1 "diag_manager/include/fms_diag_fieldbuff_update.fh" 1
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

 !> @brief This code will be used by the preprocessor to generate an implementation
 !! of the module procedure for the fieldbuff_update interface. The
  !! generated function is a wrapper calling 4D field/5D buffer version of the same.
FUNCTION fieldbuff_update_3d_r4(ofield_cfg, ofield_index_cfg, field_d, sample, &
  & ofb, ofc, bbounds, count_0d, num_elements, mask, weight1, missvalue,  &
  & field_num_threads, field_active_omp_level, issued_mask_ignore_warning, &
  & l_start, l_end, err_msg,  err_msg_local ) result( succeded )
    TYPE(fmsDiagOutfield_type), INTENT(in):: ofield_cfg  !< The fmsDiagOutfield_type object,
                                                         !! where "cfg" is short for configuration
    TYPE(fmsDiagOutfieldIndex_type) , INTENT(in) :: ofield_index_cfg !<The fmsDiagOutfieldIndex_type object,
                                                                !! where "cfg" is short for configuration
    REAL(r4_kind) , CONTIGUOUS, DIMENSION(:,:,:), INTENT(in), target :: field_d !< The input field data
    INTEGER, INTENT(in) :: sample  !< The index along the diurnal time axis
    REAL(r4_kind) , CONTIGUOUS, DIMENSION(:,:,:,:), INTENT(inout), target :: ofb !< Output Field Buffer
    REAL(r4_kind) , ALLOCATABLE, DIMENSION(:,:,:,:), INTENT(inout), target :: ofc !< Output Field Counter
    TYPE(fmsDiagIbounds_type), INTENT(inout) :: bbounds !< The array bounds of the ofb argument.
    REAL(r4_kind), INTENT(inout):: count_0d !< Normally the member of the buffer object of the same name.
    INTEGER, INTENT(inout) :: num_elements !< Used in counting updated buffer elements; Other functions (e.g. wrting
                                           !!field) may nprmalize output buffer elements with the same.
    LOGICAL, CONTIGUOUS, DIMENSION(:,:,:), INTENT(in), OPTIONAL, target:: mask !< The mask of the corresponding field.
    REAL(r4_kind), INTENT(in) :: weight1 !< Field data is multiplied by weight.
    REAL(r4_kind), INTENT(in) :: missvalue !< Buffer may be set to missvalue where mask is false.

    INTEGER, INTENT(inout) :: field_num_threads !< Number of OMP threads used processing the input field;
                                                   !!expected 1 if no OMP.
    INTEGER, INTENT(inout) :: field_active_omp_level !<The OMP active level for the input field; expected 0 if no OMP.

    LOGICAL, INTENT(inout) :: issued_mask_ignore_warning !< Will be set true if mask is ignored do ti missing values
                                                        !! not present for this field
    INTEGER, DIMENSION(3), INTENT(in)  :: l_start !< local start indices on 3 axes for regional output
    INTEGER, DIMENSION(3), INTENT(in)  :: l_end !< local end indices on 3 axes for regional output
    CHARACTER(len=*), INTENT(out),OPTIONAL::err_msg !< Possibly passed in by the caller, and sent to error handler
    CHARACTER(len=256), INTENT(out) :: err_msg_local !< Possibly set by bounds checker, and sent to error handler

    LOGICAL :: succeded !< True iff no errors encountered.

    !!TODO: Can the two dummy variables below be removed. The variables were introduced because of pointer
    !! bounds remapping to call the next routine. Note that mask is an optional argument, and ofc might not
    !! have been allocated. The program hangs if the approach of having pointers to the dummy arrays are not used.
    !! Also note that a sepeate logical field (in fmsDiagOutfield_type, but was used in
    !! the legacy diag manager also) is used to determine "if the mask was preent"!
    REAL(r4_kind) , DIMENSION(1),target :: ofc_dummy !> A target for ofc_ptr, in case ofc is not allocated
    LOGICAL, DIMENSION(1), target :: mask_dummy !> A target for mask_ptr, in case mask is not present

    !! For pointer bounds remapping
    REAL(r4_kind) , pointer, DIMENSION(:,:,:,:) :: field_ptr !< Pointer to the field
    REAL(r4_kind) , pointer, DIMENSION(:,:,:,:,:) :: ofb_ptr !< Pointer to the outfield buffer.
    REAL(r4_kind) , pointer, DIMENSION(:,:,:,:,:) :: ofc_ptr !< Pointer to the outfield counter.
    LOGICAL , pointer, DIMENSION(:,:,:,:) :: mask_ptr !< Pointer to the mask.

    !!Set all the pointers!
    field_ptr(1:size(field_d,1),1:size(field_d,2),1:size(field_d,3),1:1) => field_d
    ofb_ptr(1:size(ofb,1),1:size(ofb,2),1:size(ofb,3),1:1, 1:size(ofb,4)) => ofb

    !!Note that diag manager does not allocate the ofc in all situations
    if(allocated(ofc)) then
      ofc_ptr(1:size(ofc,1),1:size(ofc,2),1:size(ofc,3), 1:1, 1:size(ofc,4)) => ofc
    else
     ofc_ptr(1:1,1:1,1:1,1:1,1:1) => ofc_dummy
    endif

    IF  (PRESENT (mask)) THEN
      mask_ptr(1:size(mask,1),1:size(mask,2),1:size(mask,3),1:1) => mask
    ELSE
      mask_ptr(1:1,1:1,1:1,1:1) => mask_dummy
    ENDIF

    succeded = fieldbuff_update_r4 (ofield_cfg, ofield_index_cfg, field_ptr, sample, &
      & ofb_ptr, ofc_ptr, bbounds, count_0d, num_elements, mask_ptr, weight1, missvalue,  &
      & field_num_threads, field_active_omp_level, issued_mask_ignore_warning, &
      & l_start, l_end, err_msg,  err_msg_local)
  END FUNCTION fieldbuff_update_3d_r4


!> @brief This code will be used by the  preprocessor to generate an implementation
!! of the module procedure for the fieldbuff_update interface.
!! Updates elements of the running field output buffer (argument ofb)
!! and counter (argument ofc) based on the input field data array (argument field_d).
!! In general the formulas are :
!! A) ofb(l) = ofb(l) + (weight * field(l))**pow_value
!! B) ofc(l) = ofc(l) + weight
!! where l is a standing for some set of indices in multiple dimensions.
!! Note this function may set field object members active_omp_level and  num_threads.
!! TODO: (MDM) revisit passing in and need of field_num_threads and field_active_omp_level
   FUNCTION fieldbuff_update_r4 (ofield_cfg, ofield_index_cfg, field_d, sample, &
    & ofb, ofc, bbounds, count_0d, num_elements, mask, weight1, missvalue,  &
    & field_num_threads, field_active_omp_level, issued_mask_ignore_warning, &
    & l_start, l_end, err_msg,  err_msg_local ) result( succeded )
      TYPE(fmsDiagOutfield_type), INTENT(in) :: ofield_cfg  !< The fmsDiagOutfield_type object,
                                          !!where "cfg" is short for configuration
      TYPE(fmsDiagOutfieldIndex_type) , INTENT(in) :: ofield_index_cfg  !< The fmsDiagOutfieldIndex_type object,
      !!                                                                where "cfg" is short for configuration
      REAL(r4_kind) , CONTIGUOUS, DIMENSION(:,:,:,:), INTENT(in) :: field_d !< The input field data array.
      REAL(r4_kind) , DIMENSION(:,:,:,:,:), INTENT(inout) :: ofb !< Output Field Buffer
      REAL(r4_kind) , DIMENSION(:,:,:,:,:), INTENT(inout) :: ofc !< Output Field Counter
      TYPE(fmsDiagIbounds_type), INTENT(inout) :: bbounds  !< The array bounds of the ofb argument.
      INTEGER, INTENT(in) :: sample  !< The index along the diurnal time axis
      REAL(r4_kind), INTENT(inout):: count_0d !< Normally the member of the buffer object of the same name.
      INTEGER, INTENT(inout):: num_elements !< Used in counting updated buffer elements; Other functions (e.g. wrting
                                           !!field) may nprmalize output buffer elements with the same.
      LOGICAL, CONTIGUOUS, DIMENSION(:,:,:,:), INTENT(in), OPTIONAL:: mask !< The mask of the corresponding field.
      REAL(r4_kind), INTENT(in) :: weight1 !< Field data is multiplied by weight
      REAL(r4_kind), INTENT(in) :: missvalue !< Buffer may be set to missvalue where mask is false.

      INTEGER, INTENT(inout) :: field_num_threads !< Number of OMP threads used processing the input field;
                                                   !! expected 1 if no OMP.
      INTEGER, INTENT(inout)::field_active_omp_level !<The OMP active level for the input field; expected 0 if no OMP.

      LOGICAL, INTENT(inout) :: issued_mask_ignore_warning !< Will be set true if mask is ignored do ti missing values
                                                           !! not present for this field
      INTEGER, DIMENSION(3), INTENT(in)  :: l_start !< local start indices on 3 axes for regional output
      INTEGER, DIMENSION(3), INTENT(in)  :: l_end !< local end indices on 3 axes for regional output
      CHARACTER(len=*), INTENT(out),OPTIONAL::err_msg !< Possibly passed by the caller, and sent to error handler
      CHARACTER(len=256), INTENT(out) :: err_msg_local !< Possibly set by bounds checker, and sent to error handler

      INTEGER :: pow_value !< A copy of same variable in ofield_cfg
      CHARACTER(:), ALLOCATABLE :: output_name !< A copy of same variable in ofield_cfg
      CHARACTER(:), ALLOCATABLE :: field_name  !< A copy of same variable in ofield_cfg
      CHARACTER(:), ALLOCATABLE :: module_name !< A copy of same variable in ofield_cfg
      LOGICAL :: phys_window !< A copy of same variable in ofield_cfg
      LOGICAL :: need_compute !< A copy of same variable in ofield_cfg
      LOGICAL :: reduced_k_range !< A copy of same variable in ofield_cfg
      LOGICAL :: mask_variant !< A copy of same variable in ofield_cfg
      LOGICAL :: mask_present !< A copy of same variable in ofield_cfg
      LOGICAL :: missvalue_present !< A copy of same variable in ofield_cfg

      !< The indices copied directly from the ofield_index_cfg:
      INTEGER:: is, js, ks, ie, je, ke, hi, hj, f1, f2, f3, f4
      INTEGER:: ls, le !< start and end indices for the 4th dimension.

      INTEGER :: ksr, ker !< Loop indices used in reduced_k_range calculations
      INTEGER :: i, j, k, l, i1, j1, k1  !< Looping indices, derived from ofield_index_cfg:

      INTEGER :: numthreads
      INTEGER :: active_omp_level

      LOGICAL :: succeded !< True iff no errors encountered.
      CHARACTER(len=128):: error_string






      !!TODO: (MDM) Will the interface allow passing in is, ie?
      ls = 1
      le = SIZE(field_d, 4)

      ksr= l_start(3)
      ker= l_end(3)
      is = ofield_index_cfg%get_is()
      js = ofield_index_cfg%get_js()
      ks = ofield_index_cfg%get_ks()
      ie = ofield_index_cfg%get_ie()
      je = ofield_index_cfg%get_je()
      ke = ofield_index_cfg%get_ke()
      hi = ofield_index_cfg%get_hi()
      hj = ofield_index_cfg%get_hj()
      f1 = ofield_index_cfg%get_f1()
      f2 = ofield_index_cfg%get_f2()
      f3 = ofield_index_cfg%get_f3()
      f4 = ofield_index_cfg%get_f4()

      output_name = trim(ofield_cfg%get_output_name())
      field_name = trim(ofield_cfg%get_field_name())
      module_name = trim(ofield_cfg%get_module_name())
      pow_value = ofield_cfg%get_pow_value()
      phys_window = ofield_cfg%get_phys_window()
      reduced_k_range = ofield_cfg%get_reduced_k_range()
      need_compute = ofield_cfg%get_need_compute()
      mask_variant =  ofield_cfg%get_mask_variant()
      mask_present =  ofield_cfg%get_mask_present()
      missvalue_present = ofield_cfg%get_missvalue_present()

!$OMP CRITICAL
      field_num_threads = 1
      active_omp_level = 0




      numthreads = field_num_threads
      active_omp_level = field_active_omp_level
!$OMP END CRITICAL

      MASK_VAR_IF: IF ( mask_variant ) THEN
         IF ( need_compute ) THEN
            WRITE (error_string,'(a,"/",a)') TRIM(module_name), TRIM(output_name)
            IF ( fms_error_handler('diag_manager_mod::send_data_3d', 'module/output_field '//TRIM(error_string)//&
            & ', regional output NOT supported with mask_variant', err_msg)) THEN
               succeded = .FALSE.
               RETURN
            END IF
         END IF

         ! Should reduced_k_range data be supported with the mask_variant option   ?????
         ! If not, error message should be produced and the reduced_k_range loop below eliminated
         MASK_PR_1_IF: IF (mask_present ) THEN
            MISSVAL_PR_1_IF: IF ( missvalue_present ) THEN !!(section: mask_variant .eq. true + mask present)
               IF ( debug_diag_manager ) THEN
                  CALL bbounds%update_bounds(is-hi, ie-hi, js-hj, je-hj, ks, ke)
                  CALL fms_diag_check_out_of_bounds(ofb, bbounds, output_name, module_name, err_msg_local)
                  IF ( err_msg_local /= '' ) THEN
                     IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN
                        succeded = .FALSE.
                        RETURN
                     END IF
                  END IF
               END IF
               !!
               IF( numthreads>1 .AND. phys_window ) then
                  REDU_KR1_IF: IF ( reduced_k_range ) THEN
                     DO k= ksr, ker
                        k1= k - ksr + 1
                        DO j=js, je
                           DO i=is, ie
                              where ( mask(i-is+1+hi, j-js+1+hj, k, :) )
                                ofb(i-hi,j-hj,k1,:,sample) = addwf( ofb(i-hi,j-hj,k1,:,sample), &
                                 & field_d(i-is+1+hi, j-js+1+hj, k, :), weight1, pow_value)
                                ofc(i-hi,j-hj,k1,:,sample) = ofc(i-hi,j-hj,k1,:,sample) + weight1
                              END where
                           END DO
                        END DO
                     END DO
                  ELSE REDU_KR1_IF
                     DO k=ks, ke
                        DO j=js, je
                           DO i=is, ie
                              where ( mask(i-is+1+hi, j-js+1+hj, k, :) )
                                 ofb(i-hi,j-hj,k,:,sample) = addwf( ofb(i-hi,j-hj,k,:,sample) , &
                                 & field_d(i-is+1+hi,j-js+1+hj, k, :) , weight1, pow_value)
                                 ofc(i-hi,j-hj,k,:,sample) = ofc(i-hi,j-hj,k,:,sample) + weight1
                              END where
                           END DO
                        END DO
                     END DO
                  END IF REDU_KR1_IF
               ELSE
!$OMP CRITICAL
                  REDU_KR2_IF: IF ( reduced_k_range ) THEN
                     DO k= ksr, ker
                        k1= k - ksr + 1
                        DO j=js, je
                           DO i=is, ie
                              where ( mask(i-is+1+hi, j-js+1+hj, k, :) )
                                 ofb(i-hi,j-hj,k1,:,sample) = addwf( ofb(i-hi,j-hj,k1,:,sample) , &
                                 & field_d(i-is+1+hi, j-js+1+hj, k, :) , weight1, pow_value)
                                 ofc(i-hi,j-hj,k1,:,sample) = ofc(i-hi,j-hj,k1,:,sample) + weight1
                              END where
                           END DO
                        END DO
                     END DO
                  ELSE REDU_KR2_IF
                     DO k=ks, ke
                        DO j=js, je
                           DO i=is, ie
                              where ( mask(i-is+1+hi, j-js+1+hj, k, :) )
                                 ofb(i-hi,j-hj,k,:,sample) = addwf( ofb(i-hi,j-hj,k,:,sample) , &
                                 & field_d(i-is+1+hi,j-js+1+hj, k, :) , weight1, pow_value)
                                 ofc(i-hi,j-hj,k,:,sample) = ofc(i-hi,j-hj,k,:,sample) + weight1
                              END where
                           END DO
                        END DO
                     END DO
                  END IF REDU_KR2_IF
!$OMP END CRITICAL
               END IF
            ELSE MISSVAL_PR_1_IF
               WRITE (error_string,'(a,"/",a)') TRIM(module_name), TRIM(output_name)
               IF(fms_error_handler('diag_manager_mod::send_data_3d', &
               & 'module/output_field '//TRIM(error_string)//', variable mask but no missing value defined', &
               & err_msg)) THEN
                  succeded = .FALSE.
                  RETURN
               END IF
            END IF  MISSVAL_PR_1_IF
         ELSE MASK_PR_1_IF ! no mask present
            WRITE (error_string,'(a,"/",a)') TRIM(module_name), TRIM(output_name)
            IF(fms_error_handler('diag_manager_mod::send_data_3d','module/output_field '//TRIM(error_string)//&
            & ', variable mask but no mask given', err_msg)) THEN
               succeded = .FALSE.
               RETURN
            END IF
         END IF MASK_PR_1_IF
      ELSE MASK_VAR_IF
         MASK_PR_2_IF: IF (mask_present ) THEN
            MISSVAL_PR_2_IF: IF ( missvalue_present ) THEN !!section:(mask_var false +mask present +missval prsnt)
               NDCMP_RKR_1_IF: IF ( need_compute ) THEN
                  IF (numthreads>1 .AND. phys_window) then
                     DO k = l_start(3), l_end(3)
                        k1 = k-l_start(3)+1
                        DO j = js, je
                           DO i = is, ie
                              IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                              & j <= l_end(2)+hj ) THEN
                                 i1 = i-l_start(1)-hi+1
                                 j1=  j-l_start(2)-hj+1
                                 where ( mask(i-is+1+hi, j-js+1+hj, k, :) )
                                    ofb(i1,j1,k1,:,sample) = addwf( ofb(i1,j1,k1,:,sample) , &
                                    & field_d(i-is+1+hi,j-js+1+hj, k, :) , weight1, pow_value)
                                 elsewhere
                                    ofb(i1,j1,k1,:,sample) = missvalue
                                 END where
                              END IF
                           END DO
                        END DO
                     END DO
                  ELSE
!$OMP CRITICAL
                     DO k = l_start(3), l_end(3)
                        k1 = k-l_start(3)+1
                        DO j = js, je
                           DO i = is, ie
                              IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                              & j <= l_end(2)+hj ) THEN
                                 i1 = i-l_start(1)-hi+1
                                 j1=  j-l_start(2)-hj+1
                                 where ( mask(i-is+1+hi, j-js+1+hj, k, :) )
                                    ofb(i1,j1,k1,:,sample) = addwf( ofb(i1,j1,k1,:,sample) , &
                                    & field_d(i-is+1+hi,j-js+1+hj, k, :) , weight1, pow_value)
                                 elsewhere
                                    ofb(i1,j1,k1,:,sample) = missvalue
                                 END where
                              END IF
                           END DO
                        END DO
                     END DO
!$OMP END CRITICAL
                  ENDIF
!$OMP CRITICAL
                  DO l = ls, le
                    DO j = js, je
                      DO i = is, ie
                          IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                            & j <= l_end(2)+hj ) THEN
                            num_elements = num_elements + l_end(3) - l_start(3) + 1
                          END IF
                      END DO
                   END DO
                  END DO
!$OMP END CRITICAL
               ELSE IF ( reduced_k_range ) THEN NDCMP_RKR_1_IF
                  IF (numthreads>1 .AND. phys_window) then
                     DO k=ksr, ker
                        k1 = k - ksr + 1
                        DO j=js, je
                           DO i=is, ie
                              where ( mask(i-is+1+hi,j-js+1+hj, k, :) )
                                 ofb(i-hi,j-hj,k1,:,sample) = addwf( ofb(i-hi,j-hj,k1,:,sample) , &
                                 & field_d(i-is+1+hi,j-js+1+hj, k, :) , weight1, pow_value)
                              elsewhere
                                 ofb(i-hi,j-hj,k1,:,sample)= missvalue
                              END where
                           END DO
                        END DO
                     END DO
                  ELSE
!$OMP CRITICAL
                     DO k=ksr, ker
                        k1 = k - ksr + 1
                        DO j=js, je
                           DO i=is, ie
                              where ( mask(i-is+1+hi,j-js+1+hj, k, :) )
                                 ofb(i-hi,j-hj,k1,:,sample) = addwf( ofb(i-hi,j-hj,k1,:,sample) , &
                                 & field_d(i-is+1+hi,j-js+1+hj, k, :) , weight1, pow_value)
                              elsewhere
                                 ofb(i-hi,j-hj,k1,:,sample)= missvalue
                              END where
                           END DO
                        END DO
                     END DO
!$OMP END CRITICAL
                  END IF
               ELSE NDCMP_RKR_1_IF
                  IF ( debug_diag_manager ) THEN
                     CALL bbounds%update_bounds(is-hi, ie-hi, js-hj, je-hj, ks, ke)
                     CALL fms_diag_check_out_of_bounds(ofb, bbounds, output_name, module_name, err_msg_local)
                     IF ( err_msg_local /= '' ) THEN
                        IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN
                           succeded = .FALSE.
                           RETURN
                        END IF
                     END IF
                  END IF
                  IF (numthreads>1 .AND. phys_window) then
                     DO k=ks, ke
                        DO j=js, je
                           DO i=is, ie
                              where ( mask(i-is+1+hi,j-js+1+hj, k, :) )
                                 ofb(i-hi,j-hj,k,:,sample) = addwf( ofb(i-hi,j-hj,k,:,sample) , &
                                 & field_d(i-is+1+hi,j-js+1+hj, k, :) , weight1, pow_value)
                              elsewhere
                                 ofb(i-hi,j-hj,k,:,sample)= missvalue
                              END where
                           END DO
                        END DO
                     END DO
                  ELSE
!$OMP CRITICAL
                     DO k=ks, ke
                        DO j=js, je
                           DO i=is, ie
                              where ( mask(i-is+1+hi,j-js+1+hj, k, :) )
                                 ofb(i-hi,j-hj,k,:,sample) = addwf( ofb(i-hi,j-hj,k,:,sample) , &
                                 & field_d(i-is+1+hi,j-js+1+hj, k, :) , weight1, pow_value)
                              elsewhere
                                 ofb(i-hi,j-hj,k,:,sample)= missvalue
                              END where
                           END DO
                        END DO
                     END DO
!$OMP END CRITICAL
                  END IF
               END IF NDCMP_RKR_1_IF
!$OMP CRITICAL
               IF ( need_compute .AND. .NOT.phys_window ) THEN
                  IF ( ANY(mask(l_start(1)+hi:l_end(1)+hi,l_start(2)+hj:l_end(2)+hj,l_start(3):l_end(3), :)) ) &
                  count_0d = count_0d + weight1
               ELSE
                  IF ( ANY(mask(f1:f2,f3:f4,ks:ke,:)) ) count_0d = count_0d + weight1
               END IF
!$OMP END CRITICAL
            ELSE MISSVAL_PR_2_IF !! (section: mask_varian .eq. false + mask present + miss value not present)
               IF (   (.NOT.ALL(mask(f1:f2,f3:f4,ks:ke,:)) .AND. mpp_pe() .EQ. mpp_root_pe()).AND.&
               &  .NOT. issued_mask_ignore_warning) THEN
                  ! <ERROR STATUS="WARNING">
                  !   Mask will be ignored since missing values were not specified for field <field_name>
                  !   in module <module_name>
                  ! </ERROR>
                  CALL error_mesg('diag_manager_mod::send_data_3d',&
                    & 'Mask will be ignored since missing values were not specified for field '//&
                    & trim(field_name)//' in module '//&
                    & trim(module_name), WARNING)
                  issued_mask_ignore_warning = .TRUE.
               END IF
               NDCMP_RKR_2_IF: IF ( need_compute ) THEN
                  IF (numthreads>1 .AND. phys_window) then
                     DO j = js, je
                        DO i = is, ie
                           IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                           & j <= l_end(2)+hj ) THEN
                              i1 = i-l_start(1)-hi+1
                              j1 =  j-l_start(2)-hj+1
                              ofb(i1,j1,:,:,sample)= addwf(ofb(i1,j1,:,:,sample) , &
                              & field_d(i-is+1+hi,j-js+1+hj,l_start(3):l_end(3), :) , weight1, pow_value)
                           END IF
                        END DO
                     END DO
                  ELSE
!$OMP CRITICAL
                     DO j = js, je
                        DO i = is, ie
                           IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                           & j <= l_end(2)+hj ) THEN
                              i1 = i-l_start(1)-hi+1
                              j1 =  j-l_start(2)-hj+1
                              ofb(i1,j1,:,:,sample) = addwf( ofb(i1,j1,:,:,sample) , &
                              & field_d(i-is+1+hi,j-js+1+hj,l_start(3):l_end(3), :) , weight1, pow_value)
                           END IF
                        END DO
                     END DO
!$OMP END CRITICAL
                  END IF
!$OMP CRITICAL
                  DO l = ls, le
                    DO j = js, je
                      DO i = is, ie
                        IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                          & j <= l_end(2)+hj ) THEN
                            num_elements = num_elements + l_end(3)-l_start(3)+1
                        END IF
                      END DO
                    END DO
                  END DO
!$OMP END CRITICAL
               ELSE IF ( reduced_k_range ) THEN NDCMP_RKR_2_IF
                  IF (numthreads>1 .AND. phys_window) then
                     ksr= l_start(3)
                     ker= l_end(3)
                     ofb(is-hi:ie-hi,js-hj:je-hj,:,:,sample) = addwf( ofb(is-hi:ie-hi,js-hj:je-hj,:,:,sample) , &
                     & field_d(f1:f2,f3:f4,ksr:ker, :) , weight1, pow_value)
                  ELSE
!$OMP CRITICAL
                     ksr= l_start(3)
                     ker= l_end(3)
                     ofb(is-hi:ie-hi,js-hj:je-hj,:,:,sample) = addwf( ofb(is-hi:ie-hi,js-hj:je-hj,:,:,sample) , &
                     & field_d(f1:f2,f3:f4,ksr:ker,:) , weight1, pow_value)
!$OMP END CRITICAL
                  END IF
               ELSE NDCMP_RKR_2_IF
                  IF ( debug_diag_manager ) THEN
                     CALL bbounds%update_bounds(is-hi, ie-hi, js-hj, je-hj, ks, ke)
                     CALL fms_diag_check_out_of_bounds(ofb, bbounds, output_name, module_name, err_msg_local)
                     IF ( err_msg_local /= '') THEN
                        IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN
                           succeded = .FALSE.
                           RETURN
                        END IF
                     END IF
                  END IF
                  IF (numthreads>1 .AND. phys_window) then
                     ofb(is-hi:ie-hi,js-hj:je-hj,ks:ke,:,sample) =&
                     & addwf( ofb(is-hi:ie-hi,js-hj:je-hj,ks:ke,:,sample) , &
                     & field_d(f1:f2,f3:f4,ks:ke,:) , weight1, pow_value)
                  ELSE
!$OMP CRITICAL
                     ofb(is-hi:ie-hi,js-hj:je-hj,ks:ke,:,sample) =&
                     & addwf(ofb(is-hi:ie-hi,js-hj:je-hj,ks:ke,:,sample) , &
                     & field_d(f1:f2,f3:f4,ks:ke,:) , weight1, pow_value)
!$OMP END CRITICAL
                  END IF
               END IF NDCMP_RKR_2_IF
!$OMP CRITICAL
               IF ( .NOT.phys_window ) count_0d = count_0d + weight1
!$OMP END CRITICAL
            END IF MISSVAL_PR_2_IF
         ELSE MASK_PR_2_IF !!(section: mask_variant .eq. false + mask not present + missvalue)
            MISSVAL_PR_3_IF: IF (missvalue_present ) THEN
               NDCMP_RKR_3_IF: IF ( need_compute ) THEN
                  NTAPW_IF: If( numthreads>1 .AND. phys_window ) then
                     DO k = l_start(3), l_end(3)
                        k1 = k - l_start(3) + 1
                        DO j = js, je
                           DO i = is, ie
                              IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                              & j <= l_end(2)+hj) THEN
                                 i1 = i-l_start(1)-hi+1
                                 j1=  j-l_start(2)-hj+1
                                 where ( field_d(i-is+1+hi,j-js+1+hj, k, :) /= missvalue )
                                    ofb(i1,j1,k1,:,sample) = addwf( ofb(i1,j1,k1,:,sample) , &
                                    & field_d(i-is+1+hi,j-js+1+hj, k, :) , weight1, pow_value)
                                  elsewhere
                                    ofb(i1,j1,k1,:,sample) = missvalue
                                 END where
                              END IF
                           END DO
                        END DO
                     END DO
                  ELSE NTAPW_IF
!$OMP CRITICAL
                     DO k = l_start(3), l_end(3)
                        k1 = k - l_start(3) + 1
                        DO j = js, je
                           DO i = is, ie
                              IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                              & j <= l_end(2)+hj) THEN
                                 i1 = i-l_start(1)-hi+1
                                 j1=  j-l_start(2)-hj+1
                                 where ( field_d(i-is+1+hi,j-js+1+hj, k, :) /= missvalue )
                                    ofb(i1,j1,k1,:,sample) = addwf( ofb(i1,j1,k1,:,sample) , &
                                    & field_d(i-is+1+hi,j-js+1+hj, k, :) , weight1, pow_value)
                                 elsewhere
                                    ofb(i1,j1,k1,:,sample) = missvalue
                                 END where
                              END IF
                           END DO
                        END DO
                     END DO
!$OMP END CRITICAL
                  END IF NTAPW_IF
!$OMP CRITICAL
                  DO l = ls, le
                    DO j = js, je
                      DO i = is, ie
                        IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                          & j <= l_end(2)+hj) THEN
                          num_elements = num_elements + l_end(3) - l_start(3) + 1
                        END IF
                      END DO
                    END DO
                  END DO
                  IF ( .NOT.phys_window ) THEN
                     outer0: DO l = ls, le
                       DO k = l_start(3), l_end(3)
                          DO j=l_start(2)+hj, l_end(2)+hj
                             DO i=l_start(1)+hi, l_end(1)+hi
                                IF (field_d(i,j, k, l) /= missvalue ) THEN
                                  count_0d = count_0d + weight1
                                  EXIT outer0
                              END IF
                            END DO
                         END DO
                      END DO
                    END DO outer0
                  END IF
!$OMP END CRITICAL
               ELSE IF ( reduced_k_range ) THEN NDCMP_RKR_3_IF
                  if( numthreads>1 .AND. phys_window ) then
                     ksr= l_start(3)
                     ker= l_end(3)
                     DO k = ksr, ker
                        k1 = k - ksr + 1
                        DO j=js, je
                           DO i=is, ie
                              where ( field_d(i-is+1+hi,j-js+1+hj, k, :) /= missvalue )
                                 ofb(i-hi,j-hj,k1,:,sample) =  addwf(ofb(i-hi,j-hj,k1,:,sample) , &
                                 & field_d(i-is+1+hi,j-js+1+hj, k, :) , weight1, pow_value)
                              elsewhere
                                 ofb(i-hi,j-hj,k1,:,sample) = missvalue
                              END where
                           END DO
                        END DO
                     END DO
                  else
!$OMP CRITICAL
                     ksr= l_start(3)
                     ker= l_end(3)
                     DO k = ksr, ker
                        k1 = k - ksr + 1
                        DO j=js, je
                           DO i=is, ie
                              where ( field_d(i-is+1+hi,j-js+1+hj, k, :) /= missvalue )
                                 ofb(i-hi,j-hj,k1,:,sample) = addwf( ofb(i-hi,j-hj,k1,:,sample) , &
                                 & field_d(i-is+1+hi,j-js+1+hj, k, :) , weight1, pow_value)
                              elsewhere
                                 ofb(i-hi,j-hj,k1,:,sample) = missvalue
                              END where
                           END DO
                        END DO
                     END DO
!$OMP END CRITICAL
                  END IF
!$OMP CRITICAL
                  outer3: DO l = ls, le
                    DO k = ksr, ker
                      k1=k-ksr+1
                      DO j=f3, f4
                        DO i=f1, f2
                           IF ( field_d(i,j, k, l) /= missvalue ) THEN
                              count_0d = count_0d + weight1
                              EXIT outer3
                           END IF
                        END DO
                     END DO
                  END DO
                END DO outer3
!$OMP END CRITICAL
               ELSE NDCMP_RKR_3_IF
                  IF ( debug_diag_manager ) THEN
                     CALL bbounds%update_bounds(is-hi, ie-hi, js-hj, je-hj, ks, ke)
                     CALL fms_diag_check_out_of_bounds(ofb, bbounds, output_name, module_name, err_msg_local)
                     IF ( err_msg_local /= '' ) THEN
                        IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN
                           succeded = .FALSE.
                           RETURN
                        END IF
                     END IF
                  END IF
                  IF( numthreads > 1 .AND. phys_window ) then
                     DO k=ks, ke
                        DO j=js, je
                           DO i=is, ie
                              where ( field_d(i-is+1+hi,j-js+1+hj, k, :) /= missvalue )
                                 ofb(i-hi,j-hj,k,:,sample) = addwf( ofb(i-hi,j-hj,k,:,sample) , &
                                 & field_d(i-is+1+hi,j-js+1+hj, k, :) , weight1, pow_value)
                              elsewhere
                                 ofb(i-hi,j-hj,k,:,sample) = missvalue
                              END where
                           END DO
                        END DO
                     END DO
                  ELSE
!$OMP CRITICAL
                     DO k=ks, ke
                        DO j=js, je
                           DO i=is, ie
                              where ( field_d(i-is+1+hi,j-js+1+hj, k, :) /= missvalue )
                                 ofb(i-hi,j-hj,k,:,sample) = addwf( ofb(i-hi,j-hj,k,:,sample) , &
                                 & field_d(i-is+1+hi,j-js+1+hj, k, :) , weight1, pow_value)
                              elsewhere
                                 ofb(i-hi,j-hj,k,:,sample) = missvalue
                              END where
                           END DO
                        END DO
                     END DO
!$OMP END CRITICAL
                  END IF
!$OMP CRITICAL
                  outer1: DO l = ls, le
                    DO k=ks, ke
                       DO j=f3, f4
                         DO i=f1, f2
                            IF ( field_d(i,j, k, l) /= missvalue ) THEN
                              count_0d = count_0d + weight1
                              EXIT outer1
                            END IF
                          END DO
                      END DO
                    END DO
                  END DO outer1
!$OMP END CRITICAL
               END IF NDCMP_RKR_3_IF
            ELSE MISSVAL_PR_3_IF !!(section: mask_variant .eq. false + mask not present + missvalue not present)
               NDCMP_RKR_4_IF: IF ( need_compute ) THEN
                  IF( numthreads > 1 .AND. phys_window ) then
                     DO j = js, je
                        DO i = is, ie
                           IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                           & j <= l_end(2)+hj ) THEN
                              i1 = i-l_start(1)-hi+1
                              j1=  j-l_start(2)-hj+1
                              ofb(i1,j1,:,:,sample) = addwf( ofb(i1,j1,:,:,sample) , &
                              & field_d(i-is+1+hi,j-js+1+hj,l_start(3):l_end(3), :) , weight1, pow_value)
                           END IF
                        END DO
                     END DO
                  ELSE
!$OMP CRITICAL
                     DO j = js, je
                        DO i = is, ie
                           IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                           & j <= l_end(2)+hj ) THEN
                              i1 = i-l_start(1)-hi+1
                              j1=  j-l_start(2)-hj+1
                              ofb(i1,j1,:,:,sample) = addwf(ofb(i1,j1,:,:,sample), &
                                &  field_d(i-is+1+hi,j-js+1+hj,l_start(3):l_end(3), :), weight1, pow_value)
                           END IF
                        END DO
                     END DO
!$OMP END CRITICAL
                  END IF
!$OMP CRITICAL
                  DO l = ls, le
                    DO j = js, je
                      DO i = is, ie
                        IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                           & j <= l_end(2)+hj ) THEN
                          num_elements = num_elements + l_end(3)-l_start(3)+1
                        END IF
                      END DO
                    END DO
                  END DO
!$OMP END CRITICAL
                  ! Accumulate time average
               ELSE IF ( reduced_k_range ) THEN NDCMP_RKR_4_IF
                  ksr= l_start(3)
                  ker= l_end(3)
                  IF( numthreads > 1 .AND. phys_window ) then
                     ofb(is-hi:ie-hi,js-hj:je-hj,:,:,sample) =&
                     & addwf( ofb(is-hi:ie-hi,js-hj:je-hj,:,:,sample) , &
                     & field_d(f1:f2,f3:f4,ksr:ker, :) , weight1, pow_value)
                  ELSE
!$OMP CRITICAL
                     ofb(is-hi:ie-hi,js-hj:je-hj,:,:,sample) =&
                     & addwf( ofb(is-hi:ie-hi,js-hj:je-hj,:,:,sample) , &
                     & field_d(f1:f2,f3:f4,ksr:ker, :) , weight1, pow_value)
!$OMP END CRITICAL
                  END IF

               ELSE NDCMP_RKR_4_IF
                  IF ( debug_diag_manager ) THEN
                     CALL bbounds%update_bounds(is-hi, ie-hi, js-hj, je-hj, ks, ke)
                     CALL fms_diag_check_out_of_bounds(ofb, bbounds, output_name, module_name, err_msg_local)
                     IF ( err_msg_local /= '' ) THEN
                        IF (fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN
                           succeded = .FALSE.
                           RETURN
                        END IF
                     END IF
                  END IF
                  IF( numthreads > 1 .AND. phys_window ) then
                     ofb(is-hi:ie-hi,js-hj:je-hj,ks:ke,:,sample) =&
                     & addwf( ofb(is-hi:ie-hi,js-hj:je-hj,ks:ke,:,sample) , &
                     & field_d(f1:f2,f3:f4,ks:ke, :) , weight1, pow_value)
                  ELSE
!$OMP CRITICAL
                     ofb(is-hi:ie-hi,js-hj:je-hj,ks:ke,:,sample) =&
                     & addwf( ofb(is-hi:ie-hi,js-hj:je-hj,ks:ke,:,sample) , &
                     & field_d(f1:f2,f3:f4,ks:ke, :) , weight1, pow_value)
                     !!
!$OMP END CRITICAL
                  END IF
               END IF NDCMP_RKR_4_IF
!$OMP CRITICAL
               IF ( .NOT.phys_window ) count_0d = count_0d + weight1
!$OMP END CRITICAL
            END IF MISSVAL_PR_3_IF
         END IF MASK_PR_2_IF ! if mask present
      END IF MASK_VAR_IF

!$OMP CRITICAL
      IF ( .NOT.need_compute .AND. .NOT.reduced_k_range ) num_elements = num_elements + &
        & (ie-is+1)*(je-js+1)*(ke-ks+1)*(le-ls+1)
      IF ( reduced_k_range ) num_elements = num_elements + &
        & (ie-is+1)*(je-js+1)*(ker-ksr+1)*(le-ls+1)
!$OMP END CRITICAL

      succeded = .TRUE.
      RETURN

    END FUNCTION fieldbuff_update_r4


  !> @brief This code will be used by the  preprocessor to generate an implementation
  !! of the module procedure for the fieldbuff_copy_fieldvals interface. The
  !! generated function is a wrapper calling 4D field/5D buffer version of the same.
    FUNCTION fieldbuff_copy_fieldvals_3d_r4 (ofield_cfg, ofield_index_cfg, field, sample, ofb, &
      & bbounds, count_0d, mask,  missvalue,  &
      & l_start, l_end, err_msg,  err_msg_local) result( succeded )
         TYPE(fmsDiagOutfield_type), INTENT(in) :: ofield_cfg   !< The fmsDiagOutfield_type object,
                                                                !!where "cfg" is short for configuration
         TYPE(fmsDiagOutfieldIndex_type) , INTENT(in) :: ofield_index_cfg    !< The fmsDiagOutfieldIndex_type object,
                                                                          !! where "cfg" is short for configuration
         REAL(r4_kind) , CONTIGUOUS,  DIMENSION(:,:,:),INTENT(in),target:: field !< The field value array.
         INTEGER, INTENT(in) :: sample   !< index along the diurnal time axis
         REAL(r4_kind) , ALLOCATABLE, DIMENSION(:,:,:,:),INTENT(inout),target::ofb !<The Output Field Buffer
         TYPE(fmsDiagIbounds_type), INTENT(inout) :: bbounds  !< The array bounds of the ofb argument.
         REAL(r4_kind) , INTENT(inout) :: count_0d   !< Normally the member of the buffer of same name,
         LOGICAL, CONTIGUOUS, DIMENSION(:,:,:),INTENT(in),OPTIONAL,target::mask !<The mask of the corresponding field.
         REAL(r4_kind) , INTENT(in) :: missvalue !< buffer may be set to this value where mask is false.
         INTEGER, DIMENSION(3), INTENT(in)  :: l_start   !< local start indices on spatial axes for regional output
         INTEGER, DIMENSION(3), INTENT(in)  :: l_end   !< local end indices on spatial for regional output
         CHARACTER(len=*), INTENT(out),OPTIONAL::err_msg !< Possibly passed in by the caller,and sent to handler
         CHARACTER(len=256), INTENT(out) :: err_msg_local !< Possibly set by bounds checker, and sent to handler

         LOGICAL :: succeded !< True iff no errors encountered.

        !!TODO: Can the dummy variable array below be removed. The variable was introduced because of pointer
        !! bounds remapping to call the next routine. Note that mask is an optional argument.
        !! The program hangs if the approach of having a pointer to the dummy array is not used.
        !! Also note that a separate logical field (in fmsDiagOutfield_type, but was used in
        !! the legacy diag manager also) is used to determine "if the mask was preent"!
         LOGICAL, DIMENSION(1), target :: mask_dummy !> A target for mask_ptr, in case mask is not present

        !! For pointer bounds remapping:
        REAL(r4_kind) , pointer, DIMENSION(:,:,:,:) :: field_ptr!< Pointer to the field
        REAL(r4_kind) , pointer,DIMENSION(:,:,:,:,:):: ofb_ptr!< Pointer to the outfield buffer.
        LOGICAL , pointer, DIMENSION(:,:,:,:) :: mask_ptr !< Pointer to the mask.

        !Initialize all the pointers
        field_ptr(1:size(field,1),1:size(field,2),1:size(field,3),1:1) => field(:,:,:)
        ofb_ptr(1:size(ofb,1),1:size(ofb,2),1:size(ofb,3), 1:1, 1:size(ofb,4)) => ofb
        IF  (PRESENT (mask)) THEN
          mask_ptr(1:size(mask,1),1:size(mask,2),1:size(mask,3),1:1) => mask
        ELSE
          mask_ptr(1:1,1:1,1:1,1:1) => mask_dummy
        ENDIF

        succeded = fieldbuff_copy_fieldvals_r4 (ofield_cfg, ofield_index_cfg, field_ptr, sample, &
          & ofb_ptr, bbounds, count_0d, mask_ptr, missvalue,  &
          & l_start, l_end, err_msg,  err_msg_local)
    END FUNCTION fieldbuff_copy_fieldvals_3d_r4

!> @brief This code will be used by the preprocessor to generate an implementation
!! of the module procedure for the fieldbuff_copy_fieldvals interface.
!! The function may set or add to the output field buffer (argument ofb) with the input
!! field data array  (argument field)
FUNCTION fieldbuff_copy_fieldvals_r4 (ofield_cfg, ofield_index_cfg, field, sample, ofb, &
   & bbounds, count_0d, mask,  missvalue,  &
   & l_start, l_end, err_msg,  err_msg_local) result( succeded )
      TYPE(fmsDiagOutfield_type), INTENT(in) :: ofield_cfg   !< The fmsDiagOutfield_type object,
                                                             !! where "cfg" is short for configuration
      TYPE(fmsDiagOutfieldIndex_type) , INTENT(in) :: ofield_index_cfg    !< The fmsDiagOutfieldIndex_type object,
                                                               !!where "cfg" is short for configuration
      REAL(r4_kind) , CONTIGUOUS, DIMENSION(:,:,:,:), INTENT(in) :: field   !< The field value array.
      INTEGER, INTENT(in) :: sample   !< index along the diurnal time axis
      REAL(r4_kind) , DIMENSION(:,:,:,:,:), INTENT(inout) :: ofb   !< The Output Field Buffer
      TYPE(fmsDiagIbounds_type), INTENT(inout) :: bbounds !< The array bounds of the ofb argument.
      REAL(r4_kind) , INTENT(inout) :: count_0d   !< Normally the member of the buffer of same name,
      LOGICAL, CONTIGUOUS, DIMENSION(:,:,:,:), INTENT(in), OPTIONAL :: mask !< The mask of the corresponding field.
      REAL(r4_kind) , INTENT(in) :: missvalue !< buffer may be set to this value where mask is false.
      INTEGER, DIMENSION(3), INTENT(in)  :: l_start   !< local start indices on spatial axes for regional output
      INTEGER, DIMENSION(3), INTENT(in)  :: l_end   !< local end indices on spatial for regional output
      CHARACTER(len=*), INTENT(out),OPTIONAL::err_msg !< Possibly passed in by the caller, and sent to handler
      CHARACTER(len=256), INTENT(out) :: err_msg_local  !< Possibly set by bounds checker, and sent to handler
      LOGICAL :: succeded !< Return true iff errors are not encounterd.
      !!
      !!
      !< The indices copied directly from the ofield_index_cfg
      INTEGER :: is, js, ks, ie, je, ke, hi, hj, f1, f2, f3, f4

      CHARACTER(:), ALLOCATABLE  :: output_name !< A copy of same variable in ofield_cfg
      CHARACTER(:), ALLOCATABLE  :: module_name !< A copy of same variable in ofield_cfg
      LOGICAL :: need_compute !< A copy of same variable in ofield_cfg
      LOGICAL :: reduced_k_range !< A copy of same variable in ofield_cfg
      LOGICAL :: mask_present !< A copy of same variable in ofield_cfg
      LOGICAL :: missvalue_present !< A copy of same variable in ofield_cfg
      class (fmsDiagTimeReduction_type), allocatable :: time_redux !< The instance of the fmsDiagTimeReduction_type

      INTEGER :: ksr, ker !< Loop indices used in reduced_k_range calculations
      INTEGER :: i, j, k, i1, j1, k1 !< Looping indices, derived from ofield_index_cfg:
      LOGICAL :: time_max, time_min, time_sum  !< A copies of same variables in ofield_cfg%time_reduction

      ksr= l_start(3)
      ker= l_end(3)

      is = ofield_index_cfg%get_is()
      js = ofield_index_cfg%get_js()
      ks = ofield_index_cfg%get_ks()
      ie = ofield_index_cfg%get_ie()
      je = ofield_index_cfg%get_je()
      ke = ofield_index_cfg%get_ke()
      hi = ofield_index_cfg%get_hi()
      hj = ofield_index_cfg%get_hj()
      f1 = ofield_index_cfg%get_f1()
      f2 = ofield_index_cfg%get_f2()
      f3 = ofield_index_cfg%get_f3()
      f4 = ofield_index_cfg%get_f4()

      allocate(time_redux)
      call time_redux%copy(ofield_cfg%get_time_reduction())
      time_max = time_redux%is_time_max()
      time_min = time_redux%is_time_min()
      time_sum = time_redux%is_time_sum()

      output_name = trim(ofield_cfg%get_output_name())
      module_name = trim(ofield_cfg%get_module_name())
      reduced_k_range = ofield_cfg%get_reduced_k_range()
      need_compute = ofield_cfg%get_need_compute()
      mask_present =  ofield_cfg%get_mask_present()
      missvalue_present = ofield_cfg%get_missvalue_present()

         ! Add processing for Max and Min
         TIME_IF: IF ( time_max ) THEN
            MASK_PRSNT_1_IF: IF (mask_present ) THEN
               IF ( need_compute ) THEN
                  DO k = l_start(3), l_end(3)
                     k1 = k - l_start(3) + 1
                     DO j = js, je
                        DO i = is, ie
                           IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                           & j <= l_end(2)+hj ) THEN
                              i1 = i-l_start(1)-hi+1
                              j1=  j-l_start(2)-hj+1
                              WHERE ( mask(i-is+1+hi,j-js+1+hj,k,:) .AND.&
                               & field(i-is+1+hi,j-js+1+hj,k,:)>OFB(i1,j1,k1,:,sample))
                                 OFB(i1,j1,k1,:,sample) = field(i-is+1+hi,j-js+1+hj,k,:)
                              END WHERE
                           END IF
                        END DO
                     END DO
                  END DO
                  ! Maximum time value with masking
               ELSE IF ( reduced_k_range ) THEN
                  ksr = l_start(3)
                  ker = l_end(3)
                  WHERE ( mask(f1:f2,f3:f4,ksr:ker,:) .AND. &
                  & field(f1:f2,f3:f4,ksr:ker,:) > OFB(is-hi:ie-hi,js-hj:je-hj,:,:,sample))&
                  & OFB(is-hi:ie-hi,js-hj:je-hj,:,:,sample) = field(f1:f2,f3:f4,ksr:ker,:)
               ELSE
                  IF ( debug_diag_manager ) THEN
                     CALL bbounds%update_bounds(is-hi, ie-hi, js-hj, je-hj, ks, ke)
                     CALL fms_diag_check_out_of_bounds(ofb, bbounds, output_name,  module_name, err_msg_local)
                     IF ( err_msg_local /= '' ) THEN
                        IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN
                           succeded = .FALSE.
                           RETURN
                        END IF
                     END IF
                  END IF
                  WHERE ( mask(f1:f2,f3:f4,ks:ke,:) .AND.&
                  & field(f1:f2,f3:f4,ks:ke,:)>OFB(is-hi:ie-hi,js-hj:je-hj,ks:ke,:,sample))&
                  & OFB(is-hi:ie-hi,js-hj:je-hj,ks:ke,:,sample) = field(f1:f2,f3:f4,ks:ke,:)
               END IF
            ELSE MASK_PRSNT_1_IF
               IF ( need_compute ) THEN
                  DO k = l_start(3), l_end(3)
                     k1 = k - l_start(3) + 1
                     DO j = js, je
                        DO i = is, ie
                           IF(l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                           & j <= l_end(2)+hj ) THEN
                              i1 = i-l_start(1)-hi+1
                              j1 =  j-l_start(2)-hj+1
                              WHERE ( field(i-is+1+hi,j-js+1+hj,k,:) > OFB(i1,j1,k1,:,sample) )
                                 OFB(i1,j1,k1,:,sample) = field(i-is+1+hi,j-js+1+hj,k,:)
                              END WHERE
                           END IF
                        END DO
                     END DO
                  END DO
                  ! Maximum time value
               ELSE IF ( reduced_k_range ) THEN
                  ksr = l_start(3)
                  ker = l_end(3)
                  WHERE ( field(f1:f2,f3:f4,ksr:ker,:) > OFB(is-hi:ie-hi,js-hj:je-hj,:,:,sample) )&
                  & OFB(is-hi:ie-hi,js-hj:je-hj,:,:,sample) = field(f1:f2,f3:f4,ksr:ker,:)
               ELSE
                  IF ( debug_diag_manager ) THEN
                     CALL bbounds%update_bounds(is-hi, ie-hi, js-hj, je-hj, ks, ke)
                     CALL fms_diag_check_out_of_bounds(ofb, bbounds, output_name, module_name, err_msg_local)
                     IF ( err_msg_local /= '' ) THEN
                        IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN
                           succeded = .FALSE.
                           RETURN
                        END IF
                     END IF
                  END IF
                  WHERE (field(f1:f2,f3:f4,ks:ke,:) > OFB(is-hi:ie-hi,js-hj:je-hj,ks:ke,:,sample))&
                  & OFB(is-hi:ie-hi,js-hj:je-hj,ks:ke,:,sample) = field(f1:f2,f3:f4,ks:ke,:)
               END IF
            END IF MASK_PRSNT_1_IF
            count_0d = 1
            !END TIME MAX
         ELSE IF ( time_min ) THEN TiME_IF
            MASK_PRSNT_2_IF: IF (mask_present ) THEN
               IF ( need_compute ) THEN
                  DO k = l_start(3), l_end(3)
                     k1 = k - l_start(3) + 1
                     DO j = js, je
                        DO i = is, ie
                           IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                           & j <= l_end(2)+hj ) THEN
                              i1 = i-l_start(1)-hi+1
                              j1 =  j-l_start(2)-hj+1
                              WHERE ( mask(i-is+1+hi,j-js+1+hj,k,:) .AND.&
                              & field(i-is+1+hi,j-js+1+hj,k,:) < OFB(i1,j1,k1,:,sample) )
                                 OFB(i1,j1,k1,:,sample) = field(i-is+1+hi,j-js+1+hj,k,:)
                              END WHERE
                           END IF
                        END DO
                     END DO
                  END DO
                  ! Minimum time value with masking
               ELSE IF ( reduced_k_range ) THEN
                  ksr= l_start(3)
                  ker= l_end(3)
                  WHERE ( mask(f1:f2,f3:f4,ksr:ker,:) .AND.&
                  & field(f1:f2,f3:f4,ksr:ker,:) < OFB(is-hi:ie-hi,js-hj:je-hj,:,:,sample)) &
                  & OFB(is-hi:ie-hi,js-hj:je-hj,:,:,sample) = field(f1:f2,f3:f4,ksr:ker,:)
               ELSE
                  IF ( debug_diag_manager ) THEN
                     CALL bbounds%update_bounds(is-hi, ie-hi, js-hj, je-hj, ks, ke)
                     CALL fms_diag_check_out_of_bounds(ofb, bbounds, output_name, module_name, err_msg_local)
                     IF ( err_msg_local /= '' ) THEN
                        IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN
                           succeded = .FALSE.
                           RETURN
                        END IF
                     END IF
                  END IF
                  WHERE ( mask(f1:f2,f3:f4,ks:ke,:) .AND.&
                  & field(f1:f2,f3:f4,ks:ke,:) < OFB(is-hi:ie-hi,js-hj:je-hj,ks:ke,:,sample))&
                  & OFB(is-hi:ie-hi,js-hj:je-hj,ks:ke,:,sample) = field(f1:f2,f3:f4,ks:ke,:)
               END IF
            ELSE MASK_PRSNT_2_IF
               IF ( need_compute ) THEN
                  DO k = l_start(3), l_end(3)
                     k1 = k - l_start(3) + 1
                     DO j = js, je
                        DO i = is, ie
                           IF ( l_start(1)+hi <=i.AND.i<=l_end(1)+hi.AND.l_start(2)+hj<=j.AND.j<=l_end(2)+hj) THEN
                              i1 = i-l_start(1)-hi+1
                              j1=  j-l_start(2)-hj+1
                              WHERE ( field(i-is+1+hi,j-js+1+hj,k,:) < OFB(i1,j1,k1,:,sample) )
                                 OFB(i1,j1,k1,:,sample) = field(i-is+1+hi,j-js+1+hj,k,:)
                              END WHERE
                           END IF
                        END DO
                     END DO
                  END DO
                  ! Minimum time value
               ELSE IF ( reduced_k_range ) THEN
                  ksr= l_start(3)
                  ker= l_end(3)
                  WHERE ( field(f1:f2,f3:f4,ksr:ker,:) < OFB(is-hi:ie-hi,js-hj:je-hj,:,:,sample) )&
                  & OFB(is-hi:ie-hi,js-hj:je-hj,:,:,sample) = field(f1:f2,f3:f4,ksr:ker,:)
               ELSE
                  IF ( debug_diag_manager ) THEN
                     CALL bbounds%update_bounds(is-hi, ie-hi, js-hj, je-hj, ks, ke)
                     CALL fms_diag_check_out_of_bounds(ofb, bbounds, output_name, module_name, err_msg_local)
                     IF ( err_msg_local /= '' ) THEN
                        IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN
                           succeded = .FALSE.
                           RETURN
                        END IF
                     END IF
                  END IF
                  WHERE (field(f1:f2,f3:f4,ks:ke,:) < OFB(is-hi:ie-hi,js-hj:je-hj,ks:ke,:,sample))&
                  & OFB(is-hi:ie-hi,js-hj:je-hj,ks:ke,:,sample) = field(f1:f2,f3:f4,ks:ke,:)
               END IF
            END IF MASK_PRSNT_2_IF
            count_0d = 1

            !! END_TIME_MIN
         ELSE IF ( time_sum ) THEN TIME_IF
            MASK_PRSNT_3_IF: IF (mask_present ) THEN
               IF ( need_compute ) THEN
                  DO k = l_start(3), l_end(3)
                     k1 = k - l_start(3) + 1
                     DO j = js, je
                        DO i = is, ie
                           IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                           & j <= l_end(2)+hj ) THEN
                              i1 = i-l_start(1)-hi+1
                              j1 =  j-l_start(2)-hj+1
                              WHERE ( mask(i-is+1+hi,j-js+1+hj,k,:) )
                                 OFB(i1,j1,k1,:,sample) = OFB(i1,j1,k1,:,sample) + field(i-is+1+hi,j-js+1+hj,k,:)
                              END WHERE
                           END IF
                        END DO
                     END DO
                  END DO
                  ! Minimum time value with masking
               ELSE IF ( reduced_k_range ) THEN
                  ksr= l_start(3)
                  ker= l_end(3)
                  OFB(is-hi:ie-hi,js-hj:je-hj,:,:,sample) = &
                  &   OFB(is-hi:ie-hi,js-hj:je-hj,:,:,sample) + &
                  &   field(f1:f2,f3:f4,ksr:ker,:)
               ELSE
                  IF ( debug_diag_manager ) THEN
                     CALL bbounds%update_bounds(is-hi, ie-hi, js-hj, je-hj, ks, ke)
                     CALL fms_diag_check_out_of_bounds(ofb, bbounds, output_name, module_name, err_msg_local)
                     IF ( err_msg_local /= '' ) THEN
                        IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN
                           succeded = .FALSE.
                           RETURN
                        END IF
                     END IF
                  END IF
                  WHERE ( mask(f1:f2,f3:f4,ks:ke,:) ) &
                  & OFB(is-hi:ie-hi,js-hj:je-hj,ks:ke,:,sample) = &
                  &  OFB(is-hi:ie-hi,js-hj:je-hj,ks:ke,:,sample) + &
                  &  field(f1:f2,f3:f4,ks:ke,:)
               END IF
            ELSE MASK_PRSNT_3_IF
               IF ( need_compute ) THEN
                  DO k = l_start(3), l_end(3)
                     k1 = k - l_start(3) + 1
                     DO j = js, je
                        DO i = is, ie
                           IF ( l_start(1)+hi <=i.AND.i<=l_end(1)+hi.AND.l_start(2)+hj<=j.AND.j<=l_end(2)+hj) THEN
                              i1 = i-l_start(1)-hi+1
                              j1=  j-l_start(2)-hj+1
                              OFB(i1,j1,k1,:,sample) = OFB(i1,j1,k1,:,sample) + field(i-is+1+hi,j-js+1+hj,k,:)
                        END IF
                        END DO
                     END DO
                  END DO
               ELSE IF ( reduced_k_range ) THEN
                  ksr= l_start(3)
                  ker= l_end(3)
                  OFB(is-hi:ie-hi,js-hj:je-hj,:,:,sample) = OFB(is-hi:ie-hi,js-hj:je-hj,:,:,sample) + &
                  &  field(f1:f2,f3:f4,ksr:ker,:)
               ELSE
                  IF ( debug_diag_manager ) THEN
                     CALL bbounds%update_bounds(is-hi, ie-hi, js-hj, je-hj, ks, ke)
                     CALL fms_diag_check_out_of_bounds(ofb, bbounds, output_name, module_name, err_msg_local)
                     IF ( err_msg_local /= '' ) THEN
                        IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN
                           succeded = .FALSE.
                           RETURN
                        END IF
                     END IF
                  END IF
                  OFB(is-hi:ie-hi,js-hj:je-hj,ks:ke,:,sample) = OFB(is-hi:ie-hi,js-hj:je-hj,ks:ke,:,sample) + &
                  &    field(f1:f2,f3:f4,ks:ke, :)
               END IF
            END IF MASK_PRSNT_3_IF
            count_0d = 1
            !END time_sum
         ELSE TIME_IF !! ( not average, not min, not max, not sum )
            count_0d = 1
            IF ( need_compute ) THEN
               DO j = js, je
                  DO i = is, ie
                     IF (l_start(1)+hi<= i .AND. i<= l_end(1)+hi .AND. l_start(2)+hj<= j .AND. j<= l_end(2)+hj) THEN
                        i1 = i-l_start(1)-hi+1
                        j1 = j-l_start(2)-hj+1
                        OFB(i1,j1,:,:,sample) = field(i-is+1+hi,j-js+1+hj,l_start(3):l_end(3),:)
                     END IF
                  END DO
               END DO
               ! instantaneous output
            ELSE IF ( reduced_k_range ) THEN
               ksr = l_start(3)
               ker = l_end(3)
               OFB(is-hi:ie-hi,js-hj:je-hj,:,:,sample) = field(f1:f2,f3:f4,ksr:ker,:)
            ELSE
               IF ( debug_diag_manager ) THEN
                  CALL bbounds%update_bounds(is-hi, ie-hi, js-hj, je-hj, ks, ke)
                  CALL fms_diag_check_out_of_bounds(ofb, bbounds, output_name, module_name, err_msg_local)
                  IF ( err_msg_local /= '' ) THEN
                     IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN
                        succeded = .FALSE.
                        RETURN
                     END IF
                  END IF
               END IF
               OFB(is-hi:ie-hi,js-hj:je-hj,ks:ke,:,sample) = field(f1:f2,f3:f4,ks:ke,:)
            END IF

            IF (mask_present .AND. missvalue_present ) THEN
               IF ( need_compute ) THEN
                  DO k = l_start(3), l_end(3)
                     k1 = k - l_start(3) + 1
                     DO j = js, je
                        DO i = is, ie
                           IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                           & j <= l_end(2)+hj ) THEN
                              i1 = i-l_start(1)-hi+1
                              j1 =  j-l_start(2)-hj+1
                              WHERE ( .NOT.mask(i-is+1+hi,j-js+1+hj,k,:) ) &
                               & OFB(i1,j1,k1,:,sample) = missvalue
                           END IF
                        END DO
                     END DO
                  END DO
               ELSE IF ( reduced_k_range ) THEN
                  ksr= l_start(3)
                  ker= l_end(3)
                  DO k=ksr, ker
                     k1= k - ksr + 1
                     DO j=js, je
                        DO i=is, ie
                          WHERE ( mask(i-is+1+hi,j-js+1+hj,k,:) .eqv. .false.) &
                              & OFB(i-hi,j-hj,k1,:,sample)= missvalue
                        END DO
                     END DO
                  END DO
               ELSE
                  DO k=ks, ke
                     DO j=js, je
                        DO i=is, ie
                           WHERE ( .NOT. mask(i-is+1+hi,j-js+1+hj,k,:) )&
                              & OFB(i-hi,j-hj,k,:,sample)= missvalue
                        END DO
                     END DO
                  END DO
               END IF
            END IF
         END IF TIME_IF
      succeded = .TRUE.
      RETURN

   END FUNCTION fieldbuff_copy_fieldvals_r4



  !> @brief This code will be used by the  preprocessor to generate an implementation
  !! of the module procedure for the fieldbuff_copy_misvals interface. The
  !! generated function is a wrapper calling 4D field/5D buffer version of the same.
  !! TODO (MDM) the meaning of an integer rmask has to be studied.
   SUBROUTINE fieldbuff_copy_missvals_3d_r4 (ofield_cfg, ofield_index_cfg, ofb, sample, &
    & l_start, l_end, rmask, rmask_thresh, missvalue)
      TYPE(fmsDiagOutfield_type), INTENT(in) :: ofield_cfg   !< The ofield_cfg object
      TYPE(fmsDiagOutfieldIndex_type) , INTENT(in) :: ofield_index_cfg    !< The ofield_index_cfg object
      REAL(r4_kind) , CONTIGUOUS,DIMENSION(:,:,:,:),INTENT(inout),target:: ofb !<The Output Field Buffer
      INTEGER, INTENT(in) :: sample   !< index along the diurnal time axis
      REAL(r4_kind) , CONTIGUOUS,DIMENSION(:,:,:),INTENT(in),target:: rmask !<A real valued mask 4 the field
                                                      !! determined by logical mask. Updates where rmask < rmask_thresh
      REAL(r4_kind) , INTENT(in) :: rmask_thresh  !< Updates where rmask < rmask_thresh
      REAL(r4_kind) , INTENT(in) :: missvalue  !< Value used to update the buffer.
      INTEGER, DIMENSION(3), INTENT(in)  :: l_start   !< local start indices on spatial axes for regional output
      INTEGER, DIMENSION(3), INTENT(in)  :: l_end   !< local end indices on spatial for regional output

      !! These below are used in pointer bounds remapping
      REAL(r4_kind) , pointer, DIMENSION(:,:,:,:,:) :: ofb_ptr !< Pointer to the output field
                                                                        !! buffer - used in remapping
      REAL(r4_kind) , pointer, DIMENSION(:,:,:,:) :: rmask_ptr !< Pointer to the rmask - used
                                                                         !! in remapping

      !!Initialize all the pointers
      ofb_ptr(1:size(ofb,1),1:size(ofb,2),1:size(ofb,3), 1:1, 1:size(ofb,4)) => ofb
      rmask_ptr(1:size(rmask,1),1:size(rmask,2),1:size(rmask,3),1:1) => rmask

      call fieldbuff_copy_missvals_r4 (ofield_cfg, ofield_index_cfg, ofb_ptr, sample, &
        & l_start, l_end, rmask_ptr, rmask_thresh, missvalue)
  END SUBROUTINE fieldbuff_copy_missvals_3d_r4


  !> @brief This code will be used by the  preprocessor to generate an implementation
  !! of the module procedure for the fieldbuff_copy_misvals interface.
  !!  The function updates where appropriate and depending on the rmask argument,
   !! elements of the running field output buffer (argument buffer) with value missvalue.
   !! NOTE: It appears these OFB updates were introcuded by EMC MM into the tail end of the
   !! legacy send_data_3d.
   SUBROUTINE  fieldbuff_copy_missvals_r4 (ofield_cfg, ofield_index_cfg, buffer, sample, &
      & l_start, l_end, rmask, rmask_thresh, missvalue)
        TYPE(fmsDiagOutfield_type), INTENT(in) :: ofield_cfg   !< The fmsDiagOutfield_type object,
                                                               !! where "cfg" is short for configuration
        TYPE(fmsDiagOutfieldIndex_type) , INTENT(in) :: ofield_index_cfg    !< The fmsDiagOutfieldIndex_type object,
                                                               !!where "cfg" is short for configuration
        REAL(r4_kind) , INTENT(inout), DIMENSION(:,:,:,:,:) :: buffer !< the buffer to update
         INTEGER, INTENT(in) :: sample !< index along the diurnal time axis
         INTEGER, INTENT(in), DIMENSION(3):: l_start !< local start indices on 3 axes for regional output
         INTEGER, INTENT(in), DIMENSION(3):: l_end !< local end indices on 3 axes for regional output
         REAL(r4_kind) , INTENT(in), DIMENSION(:,:,:,:):: rmask  !< Updates where rmask < rmask_thresh
         REAL(r4_kind) , INTENT(in) :: rmask_thresh  !< Updates where rmask < rmask_thresh
         REAL(r4_kind) , INTENT(in) :: missvalue  !< Value used to update the buffer.

         !< Looping indices copied from corresponding one in ofield_index_cfg info:
         INTEGER :: is, js, ks, ie, je, ke, hi, hj
         !< Floags copied from corresponding one in ofield_cfg info:
         LOGICAL :: need_compute !< A copy of same variable in ofield_cfg
         LOGICAL :: reduced_k_range !< A copy of same variable in ofield_cfg
         INTEGER :: ksr, ker !< Loop indices used in reduced_k_range calculations
         !< Looping indices, derived from ofield_index_cfg info:
         INTEGER :: i, j, k, i1, j1, k1

         is = ofield_index_cfg%get_is()
         js = ofield_index_cfg%get_js()
         ks = ofield_index_cfg%get_ks()
         ie = ofield_index_cfg%get_ie()
         je = ofield_index_cfg%get_je()
         ke = ofield_index_cfg%get_ke()
         hi = ofield_index_cfg%get_hi()
         hj = ofield_index_cfg%get_hj()

         reduced_k_range = ofield_cfg%get_reduced_k_range()
         need_compute = ofield_cfg%get_need_compute()

         associate(ofb => buffer)

            ! If rmask and missing value present, then insert missing value
            IF ( need_compute ) THEN
               DO k = l_start(3), l_end(3)
                  k1 = k - l_start(3) + 1
                  DO j = js, je
                     DO i = is, ie
                        IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND.&
                        & j <= l_end(2)+hj ) THEN
                           i1 = i-l_start(1)-hi+1
                           j1 =  j-l_start(2)-hj+1
                           where ( rmask(i-is+1+hi,j-js+1+hj,k,:) <= rmask_thresh )
                              ofb(i1,j1,k1,:,sample) = missvalue
                           end where
                        END IF
                     END DO
                  END DO
               END DO
            ELSE IF ( reduced_k_range ) THEN
               ksr= l_start(3)
               ker= l_end(3)
               DO k= ksr, ker
                  k1 = k - ksr + 1
                  DO j=js, je
                     DO i=is, ie
                        where ( rmask(i-is+1+hi,j-js+1+hj,k,:) <= rmask_thresh )
                           ofb(i-hi,j-hj,k1,:,sample)= missvalue
                        endwhere
                     END DO
                  END DO
               END DO
            ELSE
               DO k=ks, ke
                  DO j=js, je
                     DO i=is, ie
                        where ( rmask(i-is+1+hi,j-js+1+hj,k,:) <= rmask_thresh )
                          ofb(i-hi,j-hj,k,:,sample)= missvalue
                        endwhere
                     END DO
                  END DO
               END DO
            END IF
         end associate
      END SUBROUTINE  fieldbuff_copy_missvals_r4
# 34 "diag_manager/include/fms_diag_fieldbuff_update.inc" 2

# 50 "diag_manager/include/fms_diag_fieldbuff_update.inc"
# 1 "diag_manager/include/fms_diag_fieldbuff_update.fh" 1
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

 !> @brief This code will be used by the preprocessor to generate an implementation
 !! of the module procedure for the fieldbuff_update interface. The
  !! generated function is a wrapper calling 4D field/5D buffer version of the same.
FUNCTION fieldbuff_update_3d_r8(ofield_cfg, ofield_index_cfg, field_d, sample, &
  & ofb, ofc, bbounds, count_0d, num_elements, mask, weight1, missvalue,  &
  & field_num_threads, field_active_omp_level, issued_mask_ignore_warning, &
  & l_start, l_end, err_msg,  err_msg_local ) result( succeded )
    TYPE(fmsDiagOutfield_type), INTENT(in):: ofield_cfg  !< The fmsDiagOutfield_type object,
                                                         !! where "cfg" is short for configuration
    TYPE(fmsDiagOutfieldIndex_type) , INTENT(in) :: ofield_index_cfg !<The fmsDiagOutfieldIndex_type object,
                                                                !! where "cfg" is short for configuration
    REAL(r8_kind) , CONTIGUOUS, DIMENSION(:,:,:), INTENT(in), target :: field_d !< The input field data
    INTEGER, INTENT(in) :: sample  !< The index along the diurnal time axis
    REAL(r8_kind) , CONTIGUOUS, DIMENSION(:,:,:,:), INTENT(inout), target :: ofb !< Output Field Buffer
    REAL(r8_kind) , ALLOCATABLE, DIMENSION(:,:,:,:), INTENT(inout), target :: ofc !< Output Field Counter
    TYPE(fmsDiagIbounds_type), INTENT(inout) :: bbounds !< The array bounds of the ofb argument.
    REAL(r8_kind), INTENT(inout):: count_0d !< Normally the member of the buffer object of the same name.
    INTEGER, INTENT(inout) :: num_elements !< Used in counting updated buffer elements; Other functions (e.g. wrting
                                           !!field) may nprmalize output buffer elements with the same.
    LOGICAL, CONTIGUOUS, DIMENSION(:,:,:), INTENT(in), OPTIONAL, target:: mask !< The mask of the corresponding field.
    REAL(r8_kind), INTENT(in) :: weight1 !< Field data is multiplied by weight.
    REAL(r8_kind), INTENT(in) :: missvalue !< Buffer may be set to missvalue where mask is false.

    INTEGER, INTENT(inout) :: field_num_threads !< Number of OMP threads used processing the input field;
                                                   !!expected 1 if no OMP.
    INTEGER, INTENT(inout) :: field_active_omp_level !<The OMP active level for the input field; expected 0 if no OMP.

    LOGICAL, INTENT(inout) :: issued_mask_ignore_warning !< Will be set true if mask is ignored do ti missing values
                                                        !! not present for this field
    INTEGER, DIMENSION(3), INTENT(in)  :: l_start !< local start indices on 3 axes for regional output
    INTEGER, DIMENSION(3), INTENT(in)  :: l_end !< local end indices on 3 axes for regional output
    CHARACTER(len=*), INTENT(out),OPTIONAL::err_msg !< Possibly passed in by the caller, and sent to error handler
    CHARACTER(len=256), INTENT(out) :: err_msg_local !< Possibly set by bounds checker, and sent to error handler

    LOGICAL :: succeded !< True iff no errors encountered.

    !!TODO: Can the two dummy variables below be removed. The variables were introduced because of pointer
    !! bounds remapping to call the next routine. Note that mask is an optional argument, and ofc might not
    !! have been allocated. The program hangs if the approach of having pointers to the dummy arrays are not used.
    !! Also note that a sepeate logical field (in fmsDiagOutfield_type, but was used in
    !! the legacy diag manager also) is used to determine "if the mask was preent"!
    REAL(r8_kind) , DIMENSION(1),target :: ofc_dummy !> A target for ofc_ptr, in case ofc is not allocated
    LOGICAL, DIMENSION(1), target :: mask_dummy !> A target for mask_ptr, in case mask is not present

    !! For pointer bounds remapping
    REAL(r8_kind) , pointer, DIMENSION(:,:,:,:) :: field_ptr !< Pointer to the field
    REAL(r8_kind) , pointer, DIMENSION(:,:,:,:,:) :: ofb_ptr !< Pointer to the outfield buffer.
    REAL(r8_kind) , pointer, DIMENSION(:,:,:,:,:) :: ofc_ptr !< Pointer to the outfield counter.
    LOGICAL , pointer, DIMENSION(:,:,:,:) :: mask_ptr !< Pointer to the mask.

    !!Set all the pointers!
    field_ptr(1:size(field_d,1),1:size(field_d,2),1:size(field_d,3),1:1) => field_d
    ofb_ptr(1:size(ofb,1),1:size(ofb,2),1:size(ofb,3),1:1, 1:size(ofb,4)) => ofb

    !!Note that diag manager does not allocate the ofc in all situations
    if(allocated(ofc)) then
      ofc_ptr(1:size(ofc,1),1:size(ofc,2),1:size(ofc,3), 1:1, 1:size(ofc,4)) => ofc
    else
     ofc_ptr(1:1,1:1,1:1,1:1,1:1) => ofc_dummy
    endif

    IF  (PRESENT (mask)) THEN
      mask_ptr(1:size(mask,1),1:size(mask,2),1:size(mask,3),1:1) => mask
    ELSE
      mask_ptr(1:1,1:1,1:1,1:1) => mask_dummy
    ENDIF

    succeded = fieldbuff_update_r8 (ofield_cfg, ofield_index_cfg, field_ptr, sample, &
      & ofb_ptr, ofc_ptr, bbounds, count_0d, num_elements, mask_ptr, weight1, missvalue,  &
      & field_num_threads, field_active_omp_level, issued_mask_ignore_warning, &
      & l_start, l_end, err_msg,  err_msg_local)
  END FUNCTION fieldbuff_update_3d_r8


!> @brief This code will be used by the  preprocessor to generate an implementation
!! of the module procedure for the fieldbuff_update interface.
!! Updates elements of the running field output buffer (argument ofb)
!! and counter (argument ofc) based on the input field data array (argument field_d).
!! In general the formulas are :
!! A) ofb(l) = ofb(l) + (weight * field(l))**pow_value
!! B) ofc(l) = ofc(l) + weight
!! where l is a standing for some set of indices in multiple dimensions.
!! Note this function may set field object members active_omp_level and  num_threads.
!! TODO: (MDM) revisit passing in and need of field_num_threads and field_active_omp_level
   FUNCTION fieldbuff_update_r8 (ofield_cfg, ofield_index_cfg, field_d, sample, &
    & ofb, ofc, bbounds, count_0d, num_elements, mask, weight1, missvalue,  &
    & field_num_threads, field_active_omp_level, issued_mask_ignore_warning, &
    & l_start, l_end, err_msg,  err_msg_local ) result( succeded )
      TYPE(fmsDiagOutfield_type), INTENT(in) :: ofield_cfg  !< The fmsDiagOutfield_type object,
                                          !!where "cfg" is short for configuration
      TYPE(fmsDiagOutfieldIndex_type) , INTENT(in) :: ofield_index_cfg  !< The fmsDiagOutfieldIndex_type object,
      !!                                                                where "cfg" is short for configuration
      REAL(r8_kind) , CONTIGUOUS, DIMENSION(:,:,:,:), INTENT(in) :: field_d !< The input field data array.
      REAL(r8_kind) , DIMENSION(:,:,:,:,:), INTENT(inout) :: ofb !< Output Field Buffer
      REAL(r8_kind) , DIMENSION(:,:,:,:,:), INTENT(inout) :: ofc !< Output Field Counter
      TYPE(fmsDiagIbounds_type), INTENT(inout) :: bbounds  !< The array bounds of the ofb argument.
      INTEGER, INTENT(in) :: sample  !< The index along the diurnal time axis
      REAL(r8_kind), INTENT(inout):: count_0d !< Normally the member of the buffer object of the same name.
      INTEGER, INTENT(inout):: num_elements !< Used in counting updated buffer elements; Other functions (e.g. wrting
                                           !!field) may nprmalize output buffer elements with the same.
      LOGICAL, CONTIGUOUS, DIMENSION(:,:,:,:), INTENT(in), OPTIONAL:: mask !< The mask of the corresponding field.
      REAL(r8_kind), INTENT(in) :: weight1 !< Field data is multiplied by weight
      REAL(r8_kind), INTENT(in) :: missvalue !< Buffer may be set to missvalue where mask is false.

      INTEGER, INTENT(inout) :: field_num_threads !< Number of OMP threads used processing the input field;
                                                   !! expected 1 if no OMP.
      INTEGER, INTENT(inout)::field_active_omp_level !<The OMP active level for the input field; expected 0 if no OMP.

      LOGICAL, INTENT(inout) :: issued_mask_ignore_warning !< Will be set true if mask is ignored do ti missing values
                                                           !! not present for this field
      INTEGER, DIMENSION(3), INTENT(in)  :: l_start !< local start indices on 3 axes for regional output
      INTEGER, DIMENSION(3), INTENT(in)  :: l_end !< local end indices on 3 axes for regional output
      CHARACTER(len=*), INTENT(out),OPTIONAL::err_msg !< Possibly passed by the caller, and sent to error handler
      CHARACTER(len=256), INTENT(out) :: err_msg_local !< Possibly set by bounds checker, and sent to error handler

      INTEGER :: pow_value !< A copy of same variable in ofield_cfg
      CHARACTER(:), ALLOCATABLE :: output_name !< A copy of same variable in ofield_cfg
      CHARACTER(:), ALLOCATABLE :: field_name  !< A copy of same variable in ofield_cfg
      CHARACTER(:), ALLOCATABLE :: module_name !< A copy of same variable in ofield_cfg
      LOGICAL :: phys_window !< A copy of same variable in ofield_cfg
      LOGICAL :: need_compute !< A copy of same variable in ofield_cfg
      LOGICAL :: reduced_k_range !< A copy of same variable in ofield_cfg
      LOGICAL :: mask_variant !< A copy of same variable in ofield_cfg
      LOGICAL :: mask_present !< A copy of same variable in ofield_cfg
      LOGICAL :: missvalue_present !< A copy of same variable in ofield_cfg

      !< The indices copied directly from the ofield_index_cfg:
      INTEGER:: is, js, ks, ie, je, ke, hi, hj, f1, f2, f3, f4
      INTEGER:: ls, le !< start and end indices for the 4th dimension.

      INTEGER :: ksr, ker !< Loop indices used in reduced_k_range calculations
      INTEGER :: i, j, k, l, i1, j1, k1  !< Looping indices, derived from ofield_index_cfg:

      INTEGER :: numthreads
      INTEGER :: active_omp_level

      LOGICAL :: succeded !< True iff no errors encountered.
      CHARACTER(len=128):: error_string






      !!TODO: (MDM) Will the interface allow passing in is, ie?
      ls = 1
      le = SIZE(field_d, 4)

      ksr= l_start(3)
      ker= l_end(3)
      is = ofield_index_cfg%get_is()
      js = ofield_index_cfg%get_js()
      ks = ofield_index_cfg%get_ks()
      ie = ofield_index_cfg%get_ie()
      je = ofield_index_cfg%get_je()
      ke = ofield_index_cfg%get_ke()
      hi = ofield_index_cfg%get_hi()
      hj = ofield_index_cfg%get_hj()
      f1 = ofield_index_cfg%get_f1()
      f2 = ofield_index_cfg%get_f2()
      f3 = ofield_index_cfg%get_f3()
      f4 = ofield_index_cfg%get_f4()

      output_name = trim(ofield_cfg%get_output_name())
      field_name = trim(ofield_cfg%get_field_name())
      module_name = trim(ofield_cfg%get_module_name())
      pow_value = ofield_cfg%get_pow_value()
      phys_window = ofield_cfg%get_phys_window()
      reduced_k_range = ofield_cfg%get_reduced_k_range()
      need_compute = ofield_cfg%get_need_compute()
      mask_variant =  ofield_cfg%get_mask_variant()
      mask_present =  ofield_cfg%get_mask_present()
      missvalue_present = ofield_cfg%get_missvalue_present()

!$OMP CRITICAL
      field_num_threads = 1
      active_omp_level = 0




      numthreads = field_num_threads
      active_omp_level = field_active_omp_level
!$OMP END CRITICAL

      MASK_VAR_IF: IF ( mask_variant ) THEN
         IF ( need_compute ) THEN
            WRITE (error_string,'(a,"/",a)') TRIM(module_name), TRIM(output_name)
            IF ( fms_error_handler('diag_manager_mod::send_data_3d', 'module/output_field '//TRIM(error_string)//&
            & ', regional output NOT supported with mask_variant', err_msg)) THEN
               succeded = .FALSE.
               RETURN
            END IF
         END IF

         ! Should reduced_k_range data be supported with the mask_variant option   ?????
         ! If not, error message should be produced and the reduced_k_range loop below eliminated
         MASK_PR_1_IF: IF (mask_present ) THEN
            MISSVAL_PR_1_IF: IF ( missvalue_present ) THEN !!(section: mask_variant .eq. true + mask present)
               IF ( debug_diag_manager ) THEN
                  CALL bbounds%update_bounds(is-hi, ie-hi, js-hj, je-hj, ks, ke)
                  CALL fms_diag_check_out_of_bounds(ofb, bbounds, output_name, module_name, err_msg_local)
                  IF ( err_msg_local /= '' ) THEN
                     IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN
                        succeded = .FALSE.
                        RETURN
                     END IF
                  END IF
               END IF
               !!
               IF( numthreads>1 .AND. phys_window ) then
                  REDU_KR1_IF: IF ( reduced_k_range ) THEN
                     DO k= ksr, ker
                        k1= k - ksr + 1
                        DO j=js, je
                           DO i=is, ie
                              where ( mask(i-is+1+hi, j-js+1+hj, k, :) )
                                ofb(i-hi,j-hj,k1,:,sample) = addwf( ofb(i-hi,j-hj,k1,:,sample), &
                                 & field_d(i-is+1+hi, j-js+1+hj, k, :), weight1, pow_value)
                                ofc(i-hi,j-hj,k1,:,sample) = ofc(i-hi,j-hj,k1,:,sample) + weight1
                              END where
                           END DO
                        END DO
                     END DO
                  ELSE REDU_KR1_IF
                     DO k=ks, ke
                        DO j=js, je
                           DO i=is, ie
                              where ( mask(i-is+1+hi, j-js+1+hj, k, :) )
                                 ofb(i-hi,j-hj,k,:,sample) = addwf( ofb(i-hi,j-hj,k,:,sample) , &
                                 & field_d(i-is+1+hi,j-js+1+hj, k, :) , weight1, pow_value)
                                 ofc(i-hi,j-hj,k,:,sample) = ofc(i-hi,j-hj,k,:,sample) + weight1
                              END where
                           END DO
                        END DO
                     END DO
                  END IF REDU_KR1_IF
               ELSE
!$OMP CRITICAL
                  REDU_KR2_IF: IF ( reduced_k_range ) THEN
                     DO k= ksr, ker
                        k1= k - ksr + 1
                        DO j=js, je
                           DO i=is, ie
                              where ( mask(i-is+1+hi, j-js+1+hj, k, :) )
                                 ofb(i-hi,j-hj,k1,:,sample) = addwf( ofb(i-hi,j-hj,k1,:,sample) , &
                                 & field_d(i-is+1+hi, j-js+1+hj, k, :) , weight1, pow_value)
                                 ofc(i-hi,j-hj,k1,:,sample) = ofc(i-hi,j-hj,k1,:,sample) + weight1
                              END where
                           END DO
                        END DO
                     END DO
                  ELSE REDU_KR2_IF
                     DO k=ks, ke
                        DO j=js, je
                           DO i=is, ie
                              where ( mask(i-is+1+hi, j-js+1+hj, k, :) )
                                 ofb(i-hi,j-hj,k,:,sample) = addwf( ofb(i-hi,j-hj,k,:,sample) , &
                                 & field_d(i-is+1+hi,j-js+1+hj, k, :) , weight1, pow_value)
                                 ofc(i-hi,j-hj,k,:,sample) = ofc(i-hi,j-hj,k,:,sample) + weight1
                              END where
                           END DO
                        END DO
                     END DO
                  END IF REDU_KR2_IF
!$OMP END CRITICAL
               END IF
            ELSE MISSVAL_PR_1_IF
               WRITE (error_string,'(a,"/",a)') TRIM(module_name), TRIM(output_name)
               IF(fms_error_handler('diag_manager_mod::send_data_3d', &
               & 'module/output_field '//TRIM(error_string)//', variable mask but no missing value defined', &
               & err_msg)) THEN
                  succeded = .FALSE.
                  RETURN
               END IF
            END IF  MISSVAL_PR_1_IF
         ELSE MASK_PR_1_IF ! no mask present
            WRITE (error_string,'(a,"/",a)') TRIM(module_name), TRIM(output_name)
            IF(fms_error_handler('diag_manager_mod::send_data_3d','module/output_field '//TRIM(error_string)//&
            & ', variable mask but no mask given', err_msg)) THEN
               succeded = .FALSE.
               RETURN
            END IF
         END IF MASK_PR_1_IF
      ELSE MASK_VAR_IF
         MASK_PR_2_IF: IF (mask_present ) THEN
            MISSVAL_PR_2_IF: IF ( missvalue_present ) THEN !!section:(mask_var false +mask present +missval prsnt)
               NDCMP_RKR_1_IF: IF ( need_compute ) THEN
                  IF (numthreads>1 .AND. phys_window) then
                     DO k = l_start(3), l_end(3)
                        k1 = k-l_start(3)+1
                        DO j = js, je
                           DO i = is, ie
                              IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                              & j <= l_end(2)+hj ) THEN
                                 i1 = i-l_start(1)-hi+1
                                 j1=  j-l_start(2)-hj+1
                                 where ( mask(i-is+1+hi, j-js+1+hj, k, :) )
                                    ofb(i1,j1,k1,:,sample) = addwf( ofb(i1,j1,k1,:,sample) , &
                                    & field_d(i-is+1+hi,j-js+1+hj, k, :) , weight1, pow_value)
                                 elsewhere
                                    ofb(i1,j1,k1,:,sample) = missvalue
                                 END where
                              END IF
                           END DO
                        END DO
                     END DO
                  ELSE
!$OMP CRITICAL
                     DO k = l_start(3), l_end(3)
                        k1 = k-l_start(3)+1
                        DO j = js, je
                           DO i = is, ie
                              IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                              & j <= l_end(2)+hj ) THEN
                                 i1 = i-l_start(1)-hi+1
                                 j1=  j-l_start(2)-hj+1
                                 where ( mask(i-is+1+hi, j-js+1+hj, k, :) )
                                    ofb(i1,j1,k1,:,sample) = addwf( ofb(i1,j1,k1,:,sample) , &
                                    & field_d(i-is+1+hi,j-js+1+hj, k, :) , weight1, pow_value)
                                 elsewhere
                                    ofb(i1,j1,k1,:,sample) = missvalue
                                 END where
                              END IF
                           END DO
                        END DO
                     END DO
!$OMP END CRITICAL
                  ENDIF
!$OMP CRITICAL
                  DO l = ls, le
                    DO j = js, je
                      DO i = is, ie
                          IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                            & j <= l_end(2)+hj ) THEN
                            num_elements = num_elements + l_end(3) - l_start(3) + 1
                          END IF
                      END DO
                   END DO
                  END DO
!$OMP END CRITICAL
               ELSE IF ( reduced_k_range ) THEN NDCMP_RKR_1_IF
                  IF (numthreads>1 .AND. phys_window) then
                     DO k=ksr, ker
                        k1 = k - ksr + 1
                        DO j=js, je
                           DO i=is, ie
                              where ( mask(i-is+1+hi,j-js+1+hj, k, :) )
                                 ofb(i-hi,j-hj,k1,:,sample) = addwf( ofb(i-hi,j-hj,k1,:,sample) , &
                                 & field_d(i-is+1+hi,j-js+1+hj, k, :) , weight1, pow_value)
                              elsewhere
                                 ofb(i-hi,j-hj,k1,:,sample)= missvalue
                              END where
                           END DO
                        END DO
                     END DO
                  ELSE
!$OMP CRITICAL
                     DO k=ksr, ker
                        k1 = k - ksr + 1
                        DO j=js, je
                           DO i=is, ie
                              where ( mask(i-is+1+hi,j-js+1+hj, k, :) )
                                 ofb(i-hi,j-hj,k1,:,sample) = addwf( ofb(i-hi,j-hj,k1,:,sample) , &
                                 & field_d(i-is+1+hi,j-js+1+hj, k, :) , weight1, pow_value)
                              elsewhere
                                 ofb(i-hi,j-hj,k1,:,sample)= missvalue
                              END where
                           END DO
                        END DO
                     END DO
!$OMP END CRITICAL
                  END IF
               ELSE NDCMP_RKR_1_IF
                  IF ( debug_diag_manager ) THEN
                     CALL bbounds%update_bounds(is-hi, ie-hi, js-hj, je-hj, ks, ke)
                     CALL fms_diag_check_out_of_bounds(ofb, bbounds, output_name, module_name, err_msg_local)
                     IF ( err_msg_local /= '' ) THEN
                        IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN
                           succeded = .FALSE.
                           RETURN
                        END IF
                     END IF
                  END IF
                  IF (numthreads>1 .AND. phys_window) then
                     DO k=ks, ke
                        DO j=js, je
                           DO i=is, ie
                              where ( mask(i-is+1+hi,j-js+1+hj, k, :) )
                                 ofb(i-hi,j-hj,k,:,sample) = addwf( ofb(i-hi,j-hj,k,:,sample) , &
                                 & field_d(i-is+1+hi,j-js+1+hj, k, :) , weight1, pow_value)
                              elsewhere
                                 ofb(i-hi,j-hj,k,:,sample)= missvalue
                              END where
                           END DO
                        END DO
                     END DO
                  ELSE
!$OMP CRITICAL
                     DO k=ks, ke
                        DO j=js, je
                           DO i=is, ie
                              where ( mask(i-is+1+hi,j-js+1+hj, k, :) )
                                 ofb(i-hi,j-hj,k,:,sample) = addwf( ofb(i-hi,j-hj,k,:,sample) , &
                                 & field_d(i-is+1+hi,j-js+1+hj, k, :) , weight1, pow_value)
                              elsewhere
                                 ofb(i-hi,j-hj,k,:,sample)= missvalue
                              END where
                           END DO
                        END DO
                     END DO
!$OMP END CRITICAL
                  END IF
               END IF NDCMP_RKR_1_IF
!$OMP CRITICAL
               IF ( need_compute .AND. .NOT.phys_window ) THEN
                  IF ( ANY(mask(l_start(1)+hi:l_end(1)+hi,l_start(2)+hj:l_end(2)+hj,l_start(3):l_end(3), :)) ) &
                  count_0d = count_0d + weight1
               ELSE
                  IF ( ANY(mask(f1:f2,f3:f4,ks:ke,:)) ) count_0d = count_0d + weight1
               END IF
!$OMP END CRITICAL
            ELSE MISSVAL_PR_2_IF !! (section: mask_varian .eq. false + mask present + miss value not present)
               IF (   (.NOT.ALL(mask(f1:f2,f3:f4,ks:ke,:)) .AND. mpp_pe() .EQ. mpp_root_pe()).AND.&
               &  .NOT. issued_mask_ignore_warning) THEN
                  ! <ERROR STATUS="WARNING">
                  !   Mask will be ignored since missing values were not specified for field <field_name>
                  !   in module <module_name>
                  ! </ERROR>
                  CALL error_mesg('diag_manager_mod::send_data_3d',&
                    & 'Mask will be ignored since missing values were not specified for field '//&
                    & trim(field_name)//' in module '//&
                    & trim(module_name), WARNING)
                  issued_mask_ignore_warning = .TRUE.
               END IF
               NDCMP_RKR_2_IF: IF ( need_compute ) THEN
                  IF (numthreads>1 .AND. phys_window) then
                     DO j = js, je
                        DO i = is, ie
                           IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                           & j <= l_end(2)+hj ) THEN
                              i1 = i-l_start(1)-hi+1
                              j1 =  j-l_start(2)-hj+1
                              ofb(i1,j1,:,:,sample)= addwf(ofb(i1,j1,:,:,sample) , &
                              & field_d(i-is+1+hi,j-js+1+hj,l_start(3):l_end(3), :) , weight1, pow_value)
                           END IF
                        END DO
                     END DO
                  ELSE
!$OMP CRITICAL
                     DO j = js, je
                        DO i = is, ie
                           IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                           & j <= l_end(2)+hj ) THEN
                              i1 = i-l_start(1)-hi+1
                              j1 =  j-l_start(2)-hj+1
                              ofb(i1,j1,:,:,sample) = addwf( ofb(i1,j1,:,:,sample) , &
                              & field_d(i-is+1+hi,j-js+1+hj,l_start(3):l_end(3), :) , weight1, pow_value)
                           END IF
                        END DO
                     END DO
!$OMP END CRITICAL
                  END IF
!$OMP CRITICAL
                  DO l = ls, le
                    DO j = js, je
                      DO i = is, ie
                        IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                          & j <= l_end(2)+hj ) THEN
                            num_elements = num_elements + l_end(3)-l_start(3)+1
                        END IF
                      END DO
                    END DO
                  END DO
!$OMP END CRITICAL
               ELSE IF ( reduced_k_range ) THEN NDCMP_RKR_2_IF
                  IF (numthreads>1 .AND. phys_window) then
                     ksr= l_start(3)
                     ker= l_end(3)
                     ofb(is-hi:ie-hi,js-hj:je-hj,:,:,sample) = addwf( ofb(is-hi:ie-hi,js-hj:je-hj,:,:,sample) , &
                     & field_d(f1:f2,f3:f4,ksr:ker, :) , weight1, pow_value)
                  ELSE
!$OMP CRITICAL
                     ksr= l_start(3)
                     ker= l_end(3)
                     ofb(is-hi:ie-hi,js-hj:je-hj,:,:,sample) = addwf( ofb(is-hi:ie-hi,js-hj:je-hj,:,:,sample) , &
                     & field_d(f1:f2,f3:f4,ksr:ker,:) , weight1, pow_value)
!$OMP END CRITICAL
                  END IF
               ELSE NDCMP_RKR_2_IF
                  IF ( debug_diag_manager ) THEN
                     CALL bbounds%update_bounds(is-hi, ie-hi, js-hj, je-hj, ks, ke)
                     CALL fms_diag_check_out_of_bounds(ofb, bbounds, output_name, module_name, err_msg_local)
                     IF ( err_msg_local /= '') THEN
                        IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN
                           succeded = .FALSE.
                           RETURN
                        END IF
                     END IF
                  END IF
                  IF (numthreads>1 .AND. phys_window) then
                     ofb(is-hi:ie-hi,js-hj:je-hj,ks:ke,:,sample) =&
                     & addwf( ofb(is-hi:ie-hi,js-hj:je-hj,ks:ke,:,sample) , &
                     & field_d(f1:f2,f3:f4,ks:ke,:) , weight1, pow_value)
                  ELSE
!$OMP CRITICAL
                     ofb(is-hi:ie-hi,js-hj:je-hj,ks:ke,:,sample) =&
                     & addwf(ofb(is-hi:ie-hi,js-hj:je-hj,ks:ke,:,sample) , &
                     & field_d(f1:f2,f3:f4,ks:ke,:) , weight1, pow_value)
!$OMP END CRITICAL
                  END IF
               END IF NDCMP_RKR_2_IF
!$OMP CRITICAL
               IF ( .NOT.phys_window ) count_0d = count_0d + weight1
!$OMP END CRITICAL
            END IF MISSVAL_PR_2_IF
         ELSE MASK_PR_2_IF !!(section: mask_variant .eq. false + mask not present + missvalue)
            MISSVAL_PR_3_IF: IF (missvalue_present ) THEN
               NDCMP_RKR_3_IF: IF ( need_compute ) THEN
                  NTAPW_IF: If( numthreads>1 .AND. phys_window ) then
                     DO k = l_start(3), l_end(3)
                        k1 = k - l_start(3) + 1
                        DO j = js, je
                           DO i = is, ie
                              IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                              & j <= l_end(2)+hj) THEN
                                 i1 = i-l_start(1)-hi+1
                                 j1=  j-l_start(2)-hj+1
                                 where ( field_d(i-is+1+hi,j-js+1+hj, k, :) /= missvalue )
                                    ofb(i1,j1,k1,:,sample) = addwf( ofb(i1,j1,k1,:,sample) , &
                                    & field_d(i-is+1+hi,j-js+1+hj, k, :) , weight1, pow_value)
                                  elsewhere
                                    ofb(i1,j1,k1,:,sample) = missvalue
                                 END where
                              END IF
                           END DO
                        END DO
                     END DO
                  ELSE NTAPW_IF
!$OMP CRITICAL
                     DO k = l_start(3), l_end(3)
                        k1 = k - l_start(3) + 1
                        DO j = js, je
                           DO i = is, ie
                              IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                              & j <= l_end(2)+hj) THEN
                                 i1 = i-l_start(1)-hi+1
                                 j1=  j-l_start(2)-hj+1
                                 where ( field_d(i-is+1+hi,j-js+1+hj, k, :) /= missvalue )
                                    ofb(i1,j1,k1,:,sample) = addwf( ofb(i1,j1,k1,:,sample) , &
                                    & field_d(i-is+1+hi,j-js+1+hj, k, :) , weight1, pow_value)
                                 elsewhere
                                    ofb(i1,j1,k1,:,sample) = missvalue
                                 END where
                              END IF
                           END DO
                        END DO
                     END DO
!$OMP END CRITICAL
                  END IF NTAPW_IF
!$OMP CRITICAL
                  DO l = ls, le
                    DO j = js, je
                      DO i = is, ie
                        IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                          & j <= l_end(2)+hj) THEN
                          num_elements = num_elements + l_end(3) - l_start(3) + 1
                        END IF
                      END DO
                    END DO
                  END DO
                  IF ( .NOT.phys_window ) THEN
                     outer0: DO l = ls, le
                       DO k = l_start(3), l_end(3)
                          DO j=l_start(2)+hj, l_end(2)+hj
                             DO i=l_start(1)+hi, l_end(1)+hi
                                IF (field_d(i,j, k, l) /= missvalue ) THEN
                                  count_0d = count_0d + weight1
                                  EXIT outer0
                              END IF
                            END DO
                         END DO
                      END DO
                    END DO outer0
                  END IF
!$OMP END CRITICAL
               ELSE IF ( reduced_k_range ) THEN NDCMP_RKR_3_IF
                  if( numthreads>1 .AND. phys_window ) then
                     ksr= l_start(3)
                     ker= l_end(3)
                     DO k = ksr, ker
                        k1 = k - ksr + 1
                        DO j=js, je
                           DO i=is, ie
                              where ( field_d(i-is+1+hi,j-js+1+hj, k, :) /= missvalue )
                                 ofb(i-hi,j-hj,k1,:,sample) =  addwf(ofb(i-hi,j-hj,k1,:,sample) , &
                                 & field_d(i-is+1+hi,j-js+1+hj, k, :) , weight1, pow_value)
                              elsewhere
                                 ofb(i-hi,j-hj,k1,:,sample) = missvalue
                              END where
                           END DO
                        END DO
                     END DO
                  else
!$OMP CRITICAL
                     ksr= l_start(3)
                     ker= l_end(3)
                     DO k = ksr, ker
                        k1 = k - ksr + 1
                        DO j=js, je
                           DO i=is, ie
                              where ( field_d(i-is+1+hi,j-js+1+hj, k, :) /= missvalue )
                                 ofb(i-hi,j-hj,k1,:,sample) = addwf( ofb(i-hi,j-hj,k1,:,sample) , &
                                 & field_d(i-is+1+hi,j-js+1+hj, k, :) , weight1, pow_value)
                              elsewhere
                                 ofb(i-hi,j-hj,k1,:,sample) = missvalue
                              END where
                           END DO
                        END DO
                     END DO
!$OMP END CRITICAL
                  END IF
!$OMP CRITICAL
                  outer3: DO l = ls, le
                    DO k = ksr, ker
                      k1=k-ksr+1
                      DO j=f3, f4
                        DO i=f1, f2
                           IF ( field_d(i,j, k, l) /= missvalue ) THEN
                              count_0d = count_0d + weight1
                              EXIT outer3
                           END IF
                        END DO
                     END DO
                  END DO
                END DO outer3
!$OMP END CRITICAL
               ELSE NDCMP_RKR_3_IF
                  IF ( debug_diag_manager ) THEN
                     CALL bbounds%update_bounds(is-hi, ie-hi, js-hj, je-hj, ks, ke)
                     CALL fms_diag_check_out_of_bounds(ofb, bbounds, output_name, module_name, err_msg_local)
                     IF ( err_msg_local /= '' ) THEN
                        IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN
                           succeded = .FALSE.
                           RETURN
                        END IF
                     END IF
                  END IF
                  IF( numthreads > 1 .AND. phys_window ) then
                     DO k=ks, ke
                        DO j=js, je
                           DO i=is, ie
                              where ( field_d(i-is+1+hi,j-js+1+hj, k, :) /= missvalue )
                                 ofb(i-hi,j-hj,k,:,sample) = addwf( ofb(i-hi,j-hj,k,:,sample) , &
                                 & field_d(i-is+1+hi,j-js+1+hj, k, :) , weight1, pow_value)
                              elsewhere
                                 ofb(i-hi,j-hj,k,:,sample) = missvalue
                              END where
                           END DO
                        END DO
                     END DO
                  ELSE
!$OMP CRITICAL
                     DO k=ks, ke
                        DO j=js, je
                           DO i=is, ie
                              where ( field_d(i-is+1+hi,j-js+1+hj, k, :) /= missvalue )
                                 ofb(i-hi,j-hj,k,:,sample) = addwf( ofb(i-hi,j-hj,k,:,sample) , &
                                 & field_d(i-is+1+hi,j-js+1+hj, k, :) , weight1, pow_value)
                              elsewhere
                                 ofb(i-hi,j-hj,k,:,sample) = missvalue
                              END where
                           END DO
                        END DO
                     END DO
!$OMP END CRITICAL
                  END IF
!$OMP CRITICAL
                  outer1: DO l = ls, le
                    DO k=ks, ke
                       DO j=f3, f4
                         DO i=f1, f2
                            IF ( field_d(i,j, k, l) /= missvalue ) THEN
                              count_0d = count_0d + weight1
                              EXIT outer1
                            END IF
                          END DO
                      END DO
                    END DO
                  END DO outer1
!$OMP END CRITICAL
               END IF NDCMP_RKR_3_IF
            ELSE MISSVAL_PR_3_IF !!(section: mask_variant .eq. false + mask not present + missvalue not present)
               NDCMP_RKR_4_IF: IF ( need_compute ) THEN
                  IF( numthreads > 1 .AND. phys_window ) then
                     DO j = js, je
                        DO i = is, ie
                           IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                           & j <= l_end(2)+hj ) THEN
                              i1 = i-l_start(1)-hi+1
                              j1=  j-l_start(2)-hj+1
                              ofb(i1,j1,:,:,sample) = addwf( ofb(i1,j1,:,:,sample) , &
                              & field_d(i-is+1+hi,j-js+1+hj,l_start(3):l_end(3), :) , weight1, pow_value)
                           END IF
                        END DO
                     END DO
                  ELSE
!$OMP CRITICAL
                     DO j = js, je
                        DO i = is, ie
                           IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                           & j <= l_end(2)+hj ) THEN
                              i1 = i-l_start(1)-hi+1
                              j1=  j-l_start(2)-hj+1
                              ofb(i1,j1,:,:,sample) = addwf(ofb(i1,j1,:,:,sample), &
                                &  field_d(i-is+1+hi,j-js+1+hj,l_start(3):l_end(3), :), weight1, pow_value)
                           END IF
                        END DO
                     END DO
!$OMP END CRITICAL
                  END IF
!$OMP CRITICAL
                  DO l = ls, le
                    DO j = js, je
                      DO i = is, ie
                        IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                           & j <= l_end(2)+hj ) THEN
                          num_elements = num_elements + l_end(3)-l_start(3)+1
                        END IF
                      END DO
                    END DO
                  END DO
!$OMP END CRITICAL
                  ! Accumulate time average
               ELSE IF ( reduced_k_range ) THEN NDCMP_RKR_4_IF
                  ksr= l_start(3)
                  ker= l_end(3)
                  IF( numthreads > 1 .AND. phys_window ) then
                     ofb(is-hi:ie-hi,js-hj:je-hj,:,:,sample) =&
                     & addwf( ofb(is-hi:ie-hi,js-hj:je-hj,:,:,sample) , &
                     & field_d(f1:f2,f3:f4,ksr:ker, :) , weight1, pow_value)
                  ELSE
!$OMP CRITICAL
                     ofb(is-hi:ie-hi,js-hj:je-hj,:,:,sample) =&
                     & addwf( ofb(is-hi:ie-hi,js-hj:je-hj,:,:,sample) , &
                     & field_d(f1:f2,f3:f4,ksr:ker, :) , weight1, pow_value)
!$OMP END CRITICAL
                  END IF

               ELSE NDCMP_RKR_4_IF
                  IF ( debug_diag_manager ) THEN
                     CALL bbounds%update_bounds(is-hi, ie-hi, js-hj, je-hj, ks, ke)
                     CALL fms_diag_check_out_of_bounds(ofb, bbounds, output_name, module_name, err_msg_local)
                     IF ( err_msg_local /= '' ) THEN
                        IF (fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN
                           succeded = .FALSE.
                           RETURN
                        END IF
                     END IF
                  END IF
                  IF( numthreads > 1 .AND. phys_window ) then
                     ofb(is-hi:ie-hi,js-hj:je-hj,ks:ke,:,sample) =&
                     & addwf( ofb(is-hi:ie-hi,js-hj:je-hj,ks:ke,:,sample) , &
                     & field_d(f1:f2,f3:f4,ks:ke, :) , weight1, pow_value)
                  ELSE
!$OMP CRITICAL
                     ofb(is-hi:ie-hi,js-hj:je-hj,ks:ke,:,sample) =&
                     & addwf( ofb(is-hi:ie-hi,js-hj:je-hj,ks:ke,:,sample) , &
                     & field_d(f1:f2,f3:f4,ks:ke, :) , weight1, pow_value)
                     !!
!$OMP END CRITICAL
                  END IF
               END IF NDCMP_RKR_4_IF
!$OMP CRITICAL
               IF ( .NOT.phys_window ) count_0d = count_0d + weight1
!$OMP END CRITICAL
            END IF MISSVAL_PR_3_IF
         END IF MASK_PR_2_IF ! if mask present
      END IF MASK_VAR_IF

!$OMP CRITICAL
      IF ( .NOT.need_compute .AND. .NOT.reduced_k_range ) num_elements = num_elements + &
        & (ie-is+1)*(je-js+1)*(ke-ks+1)*(le-ls+1)
      IF ( reduced_k_range ) num_elements = num_elements + &
        & (ie-is+1)*(je-js+1)*(ker-ksr+1)*(le-ls+1)
!$OMP END CRITICAL

      succeded = .TRUE.
      RETURN

    END FUNCTION fieldbuff_update_r8


  !> @brief This code will be used by the  preprocessor to generate an implementation
  !! of the module procedure for the fieldbuff_copy_fieldvals interface. The
  !! generated function is a wrapper calling 4D field/5D buffer version of the same.
    FUNCTION fieldbuff_copy_fieldvals_3d_r8 (ofield_cfg, ofield_index_cfg, field, sample, ofb, &
      & bbounds, count_0d, mask,  missvalue,  &
      & l_start, l_end, err_msg,  err_msg_local) result( succeded )
         TYPE(fmsDiagOutfield_type), INTENT(in) :: ofield_cfg   !< The fmsDiagOutfield_type object,
                                                                !!where "cfg" is short for configuration
         TYPE(fmsDiagOutfieldIndex_type) , INTENT(in) :: ofield_index_cfg    !< The fmsDiagOutfieldIndex_type object,
                                                                          !! where "cfg" is short for configuration
         REAL(r8_kind) , CONTIGUOUS,  DIMENSION(:,:,:),INTENT(in),target:: field !< The field value array.
         INTEGER, INTENT(in) :: sample   !< index along the diurnal time axis
         REAL(r8_kind) , ALLOCATABLE, DIMENSION(:,:,:,:),INTENT(inout),target::ofb !<The Output Field Buffer
         TYPE(fmsDiagIbounds_type), INTENT(inout) :: bbounds  !< The array bounds of the ofb argument.
         REAL(r8_kind) , INTENT(inout) :: count_0d   !< Normally the member of the buffer of same name,
         LOGICAL, CONTIGUOUS, DIMENSION(:,:,:),INTENT(in),OPTIONAL,target::mask !<The mask of the corresponding field.
         REAL(r8_kind) , INTENT(in) :: missvalue !< buffer may be set to this value where mask is false.
         INTEGER, DIMENSION(3), INTENT(in)  :: l_start   !< local start indices on spatial axes for regional output
         INTEGER, DIMENSION(3), INTENT(in)  :: l_end   !< local end indices on spatial for regional output
         CHARACTER(len=*), INTENT(out),OPTIONAL::err_msg !< Possibly passed in by the caller,and sent to handler
         CHARACTER(len=256), INTENT(out) :: err_msg_local !< Possibly set by bounds checker, and sent to handler

         LOGICAL :: succeded !< True iff no errors encountered.

        !!TODO: Can the dummy variable array below be removed. The variable was introduced because of pointer
        !! bounds remapping to call the next routine. Note that mask is an optional argument.
        !! The program hangs if the approach of having a pointer to the dummy array is not used.
        !! Also note that a separate logical field (in fmsDiagOutfield_type, but was used in
        !! the legacy diag manager also) is used to determine "if the mask was preent"!
         LOGICAL, DIMENSION(1), target :: mask_dummy !> A target for mask_ptr, in case mask is not present

        !! For pointer bounds remapping:
        REAL(r8_kind) , pointer, DIMENSION(:,:,:,:) :: field_ptr!< Pointer to the field
        REAL(r8_kind) , pointer,DIMENSION(:,:,:,:,:):: ofb_ptr!< Pointer to the outfield buffer.
        LOGICAL , pointer, DIMENSION(:,:,:,:) :: mask_ptr !< Pointer to the mask.

        !Initialize all the pointers
        field_ptr(1:size(field,1),1:size(field,2),1:size(field,3),1:1) => field(:,:,:)
        ofb_ptr(1:size(ofb,1),1:size(ofb,2),1:size(ofb,3), 1:1, 1:size(ofb,4)) => ofb
        IF  (PRESENT (mask)) THEN
          mask_ptr(1:size(mask,1),1:size(mask,2),1:size(mask,3),1:1) => mask
        ELSE
          mask_ptr(1:1,1:1,1:1,1:1) => mask_dummy
        ENDIF

        succeded = fieldbuff_copy_fieldvals_r8 (ofield_cfg, ofield_index_cfg, field_ptr, sample, &
          & ofb_ptr, bbounds, count_0d, mask_ptr, missvalue,  &
          & l_start, l_end, err_msg,  err_msg_local)
    END FUNCTION fieldbuff_copy_fieldvals_3d_r8

!> @brief This code will be used by the preprocessor to generate an implementation
!! of the module procedure for the fieldbuff_copy_fieldvals interface.
!! The function may set or add to the output field buffer (argument ofb) with the input
!! field data array  (argument field)
FUNCTION fieldbuff_copy_fieldvals_r8 (ofield_cfg, ofield_index_cfg, field, sample, ofb, &
   & bbounds, count_0d, mask,  missvalue,  &
   & l_start, l_end, err_msg,  err_msg_local) result( succeded )
      TYPE(fmsDiagOutfield_type), INTENT(in) :: ofield_cfg   !< The fmsDiagOutfield_type object,
                                                             !! where "cfg" is short for configuration
      TYPE(fmsDiagOutfieldIndex_type) , INTENT(in) :: ofield_index_cfg    !< The fmsDiagOutfieldIndex_type object,
                                                               !!where "cfg" is short for configuration
      REAL(r8_kind) , CONTIGUOUS, DIMENSION(:,:,:,:), INTENT(in) :: field   !< The field value array.
      INTEGER, INTENT(in) :: sample   !< index along the diurnal time axis
      REAL(r8_kind) , DIMENSION(:,:,:,:,:), INTENT(inout) :: ofb   !< The Output Field Buffer
      TYPE(fmsDiagIbounds_type), INTENT(inout) :: bbounds !< The array bounds of the ofb argument.
      REAL(r8_kind) , INTENT(inout) :: count_0d   !< Normally the member of the buffer of same name,
      LOGICAL, CONTIGUOUS, DIMENSION(:,:,:,:), INTENT(in), OPTIONAL :: mask !< The mask of the corresponding field.
      REAL(r8_kind) , INTENT(in) :: missvalue !< buffer may be set to this value where mask is false.
      INTEGER, DIMENSION(3), INTENT(in)  :: l_start   !< local start indices on spatial axes for regional output
      INTEGER, DIMENSION(3), INTENT(in)  :: l_end   !< local end indices on spatial for regional output
      CHARACTER(len=*), INTENT(out),OPTIONAL::err_msg !< Possibly passed in by the caller, and sent to handler
      CHARACTER(len=256), INTENT(out) :: err_msg_local  !< Possibly set by bounds checker, and sent to handler
      LOGICAL :: succeded !< Return true iff errors are not encounterd.
      !!
      !!
      !< The indices copied directly from the ofield_index_cfg
      INTEGER :: is, js, ks, ie, je, ke, hi, hj, f1, f2, f3, f4

      CHARACTER(:), ALLOCATABLE  :: output_name !< A copy of same variable in ofield_cfg
      CHARACTER(:), ALLOCATABLE  :: module_name !< A copy of same variable in ofield_cfg
      LOGICAL :: need_compute !< A copy of same variable in ofield_cfg
      LOGICAL :: reduced_k_range !< A copy of same variable in ofield_cfg
      LOGICAL :: mask_present !< A copy of same variable in ofield_cfg
      LOGICAL :: missvalue_present !< A copy of same variable in ofield_cfg
      class (fmsDiagTimeReduction_type), allocatable :: time_redux !< The instance of the fmsDiagTimeReduction_type

      INTEGER :: ksr, ker !< Loop indices used in reduced_k_range calculations
      INTEGER :: i, j, k, i1, j1, k1 !< Looping indices, derived from ofield_index_cfg:
      LOGICAL :: time_max, time_min, time_sum  !< A copies of same variables in ofield_cfg%time_reduction

      ksr= l_start(3)
      ker= l_end(3)

      is = ofield_index_cfg%get_is()
      js = ofield_index_cfg%get_js()
      ks = ofield_index_cfg%get_ks()
      ie = ofield_index_cfg%get_ie()
      je = ofield_index_cfg%get_je()
      ke = ofield_index_cfg%get_ke()
      hi = ofield_index_cfg%get_hi()
      hj = ofield_index_cfg%get_hj()
      f1 = ofield_index_cfg%get_f1()
      f2 = ofield_index_cfg%get_f2()
      f3 = ofield_index_cfg%get_f3()
      f4 = ofield_index_cfg%get_f4()

      allocate(time_redux)
      call time_redux%copy(ofield_cfg%get_time_reduction())
      time_max = time_redux%is_time_max()
      time_min = time_redux%is_time_min()
      time_sum = time_redux%is_time_sum()

      output_name = trim(ofield_cfg%get_output_name())
      module_name = trim(ofield_cfg%get_module_name())
      reduced_k_range = ofield_cfg%get_reduced_k_range()
      need_compute = ofield_cfg%get_need_compute()
      mask_present =  ofield_cfg%get_mask_present()
      missvalue_present = ofield_cfg%get_missvalue_present()

         ! Add processing for Max and Min
         TIME_IF: IF ( time_max ) THEN
            MASK_PRSNT_1_IF: IF (mask_present ) THEN
               IF ( need_compute ) THEN
                  DO k = l_start(3), l_end(3)
                     k1 = k - l_start(3) + 1
                     DO j = js, je
                        DO i = is, ie
                           IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                           & j <= l_end(2)+hj ) THEN
                              i1 = i-l_start(1)-hi+1
                              j1=  j-l_start(2)-hj+1
                              WHERE ( mask(i-is+1+hi,j-js+1+hj,k,:) .AND.&
                               & field(i-is+1+hi,j-js+1+hj,k,:)>OFB(i1,j1,k1,:,sample))
                                 OFB(i1,j1,k1,:,sample) = field(i-is+1+hi,j-js+1+hj,k,:)
                              END WHERE
                           END IF
                        END DO
                     END DO
                  END DO
                  ! Maximum time value with masking
               ELSE IF ( reduced_k_range ) THEN
                  ksr = l_start(3)
                  ker = l_end(3)
                  WHERE ( mask(f1:f2,f3:f4,ksr:ker,:) .AND. &
                  & field(f1:f2,f3:f4,ksr:ker,:) > OFB(is-hi:ie-hi,js-hj:je-hj,:,:,sample))&
                  & OFB(is-hi:ie-hi,js-hj:je-hj,:,:,sample) = field(f1:f2,f3:f4,ksr:ker,:)
               ELSE
                  IF ( debug_diag_manager ) THEN
                     CALL bbounds%update_bounds(is-hi, ie-hi, js-hj, je-hj, ks, ke)
                     CALL fms_diag_check_out_of_bounds(ofb, bbounds, output_name,  module_name, err_msg_local)
                     IF ( err_msg_local /= '' ) THEN
                        IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN
                           succeded = .FALSE.
                           RETURN
                        END IF
                     END IF
                  END IF
                  WHERE ( mask(f1:f2,f3:f4,ks:ke,:) .AND.&
                  & field(f1:f2,f3:f4,ks:ke,:)>OFB(is-hi:ie-hi,js-hj:je-hj,ks:ke,:,sample))&
                  & OFB(is-hi:ie-hi,js-hj:je-hj,ks:ke,:,sample) = field(f1:f2,f3:f4,ks:ke,:)
               END IF
            ELSE MASK_PRSNT_1_IF
               IF ( need_compute ) THEN
                  DO k = l_start(3), l_end(3)
                     k1 = k - l_start(3) + 1
                     DO j = js, je
                        DO i = is, ie
                           IF(l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                           & j <= l_end(2)+hj ) THEN
                              i1 = i-l_start(1)-hi+1
                              j1 =  j-l_start(2)-hj+1
                              WHERE ( field(i-is+1+hi,j-js+1+hj,k,:) > OFB(i1,j1,k1,:,sample) )
                                 OFB(i1,j1,k1,:,sample) = field(i-is+1+hi,j-js+1+hj,k,:)
                              END WHERE
                           END IF
                        END DO
                     END DO
                  END DO
                  ! Maximum time value
               ELSE IF ( reduced_k_range ) THEN
                  ksr = l_start(3)
                  ker = l_end(3)
                  WHERE ( field(f1:f2,f3:f4,ksr:ker,:) > OFB(is-hi:ie-hi,js-hj:je-hj,:,:,sample) )&
                  & OFB(is-hi:ie-hi,js-hj:je-hj,:,:,sample) = field(f1:f2,f3:f4,ksr:ker,:)
               ELSE
                  IF ( debug_diag_manager ) THEN
                     CALL bbounds%update_bounds(is-hi, ie-hi, js-hj, je-hj, ks, ke)
                     CALL fms_diag_check_out_of_bounds(ofb, bbounds, output_name, module_name, err_msg_local)
                     IF ( err_msg_local /= '' ) THEN
                        IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN
                           succeded = .FALSE.
                           RETURN
                        END IF
                     END IF
                  END IF
                  WHERE (field(f1:f2,f3:f4,ks:ke,:) > OFB(is-hi:ie-hi,js-hj:je-hj,ks:ke,:,sample))&
                  & OFB(is-hi:ie-hi,js-hj:je-hj,ks:ke,:,sample) = field(f1:f2,f3:f4,ks:ke,:)
               END IF
            END IF MASK_PRSNT_1_IF
            count_0d = 1
            !END TIME MAX
         ELSE IF ( time_min ) THEN TiME_IF
            MASK_PRSNT_2_IF: IF (mask_present ) THEN
               IF ( need_compute ) THEN
                  DO k = l_start(3), l_end(3)
                     k1 = k - l_start(3) + 1
                     DO j = js, je
                        DO i = is, ie
                           IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                           & j <= l_end(2)+hj ) THEN
                              i1 = i-l_start(1)-hi+1
                              j1 =  j-l_start(2)-hj+1
                              WHERE ( mask(i-is+1+hi,j-js+1+hj,k,:) .AND.&
                              & field(i-is+1+hi,j-js+1+hj,k,:) < OFB(i1,j1,k1,:,sample) )
                                 OFB(i1,j1,k1,:,sample) = field(i-is+1+hi,j-js+1+hj,k,:)
                              END WHERE
                           END IF
                        END DO
                     END DO
                  END DO
                  ! Minimum time value with masking
               ELSE IF ( reduced_k_range ) THEN
                  ksr= l_start(3)
                  ker= l_end(3)
                  WHERE ( mask(f1:f2,f3:f4,ksr:ker,:) .AND.&
                  & field(f1:f2,f3:f4,ksr:ker,:) < OFB(is-hi:ie-hi,js-hj:je-hj,:,:,sample)) &
                  & OFB(is-hi:ie-hi,js-hj:je-hj,:,:,sample) = field(f1:f2,f3:f4,ksr:ker,:)
               ELSE
                  IF ( debug_diag_manager ) THEN
                     CALL bbounds%update_bounds(is-hi, ie-hi, js-hj, je-hj, ks, ke)
                     CALL fms_diag_check_out_of_bounds(ofb, bbounds, output_name, module_name, err_msg_local)
                     IF ( err_msg_local /= '' ) THEN
                        IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN
                           succeded = .FALSE.
                           RETURN
                        END IF
                     END IF
                  END IF
                  WHERE ( mask(f1:f2,f3:f4,ks:ke,:) .AND.&
                  & field(f1:f2,f3:f4,ks:ke,:) < OFB(is-hi:ie-hi,js-hj:je-hj,ks:ke,:,sample))&
                  & OFB(is-hi:ie-hi,js-hj:je-hj,ks:ke,:,sample) = field(f1:f2,f3:f4,ks:ke,:)
               END IF
            ELSE MASK_PRSNT_2_IF
               IF ( need_compute ) THEN
                  DO k = l_start(3), l_end(3)
                     k1 = k - l_start(3) + 1
                     DO j = js, je
                        DO i = is, ie
                           IF ( l_start(1)+hi <=i.AND.i<=l_end(1)+hi.AND.l_start(2)+hj<=j.AND.j<=l_end(2)+hj) THEN
                              i1 = i-l_start(1)-hi+1
                              j1=  j-l_start(2)-hj+1
                              WHERE ( field(i-is+1+hi,j-js+1+hj,k,:) < OFB(i1,j1,k1,:,sample) )
                                 OFB(i1,j1,k1,:,sample) = field(i-is+1+hi,j-js+1+hj,k,:)
                              END WHERE
                           END IF
                        END DO
                     END DO
                  END DO
                  ! Minimum time value
               ELSE IF ( reduced_k_range ) THEN
                  ksr= l_start(3)
                  ker= l_end(3)
                  WHERE ( field(f1:f2,f3:f4,ksr:ker,:) < OFB(is-hi:ie-hi,js-hj:je-hj,:,:,sample) )&
                  & OFB(is-hi:ie-hi,js-hj:je-hj,:,:,sample) = field(f1:f2,f3:f4,ksr:ker,:)
               ELSE
                  IF ( debug_diag_manager ) THEN
                     CALL bbounds%update_bounds(is-hi, ie-hi, js-hj, je-hj, ks, ke)
                     CALL fms_diag_check_out_of_bounds(ofb, bbounds, output_name, module_name, err_msg_local)
                     IF ( err_msg_local /= '' ) THEN
                        IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN
                           succeded = .FALSE.
                           RETURN
                        END IF
                     END IF
                  END IF
                  WHERE (field(f1:f2,f3:f4,ks:ke,:) < OFB(is-hi:ie-hi,js-hj:je-hj,ks:ke,:,sample))&
                  & OFB(is-hi:ie-hi,js-hj:je-hj,ks:ke,:,sample) = field(f1:f2,f3:f4,ks:ke,:)
               END IF
            END IF MASK_PRSNT_2_IF
            count_0d = 1

            !! END_TIME_MIN
         ELSE IF ( time_sum ) THEN TIME_IF
            MASK_PRSNT_3_IF: IF (mask_present ) THEN
               IF ( need_compute ) THEN
                  DO k = l_start(3), l_end(3)
                     k1 = k - l_start(3) + 1
                     DO j = js, je
                        DO i = is, ie
                           IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                           & j <= l_end(2)+hj ) THEN
                              i1 = i-l_start(1)-hi+1
                              j1 =  j-l_start(2)-hj+1
                              WHERE ( mask(i-is+1+hi,j-js+1+hj,k,:) )
                                 OFB(i1,j1,k1,:,sample) = OFB(i1,j1,k1,:,sample) + field(i-is+1+hi,j-js+1+hj,k,:)
                              END WHERE
                           END IF
                        END DO
                     END DO
                  END DO
                  ! Minimum time value with masking
               ELSE IF ( reduced_k_range ) THEN
                  ksr= l_start(3)
                  ker= l_end(3)
                  OFB(is-hi:ie-hi,js-hj:je-hj,:,:,sample) = &
                  &   OFB(is-hi:ie-hi,js-hj:je-hj,:,:,sample) + &
                  &   field(f1:f2,f3:f4,ksr:ker,:)
               ELSE
                  IF ( debug_diag_manager ) THEN
                     CALL bbounds%update_bounds(is-hi, ie-hi, js-hj, je-hj, ks, ke)
                     CALL fms_diag_check_out_of_bounds(ofb, bbounds, output_name, module_name, err_msg_local)
                     IF ( err_msg_local /= '' ) THEN
                        IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN
                           succeded = .FALSE.
                           RETURN
                        END IF
                     END IF
                  END IF
                  WHERE ( mask(f1:f2,f3:f4,ks:ke,:) ) &
                  & OFB(is-hi:ie-hi,js-hj:je-hj,ks:ke,:,sample) = &
                  &  OFB(is-hi:ie-hi,js-hj:je-hj,ks:ke,:,sample) + &
                  &  field(f1:f2,f3:f4,ks:ke,:)
               END IF
            ELSE MASK_PRSNT_3_IF
               IF ( need_compute ) THEN
                  DO k = l_start(3), l_end(3)
                     k1 = k - l_start(3) + 1
                     DO j = js, je
                        DO i = is, ie
                           IF ( l_start(1)+hi <=i.AND.i<=l_end(1)+hi.AND.l_start(2)+hj<=j.AND.j<=l_end(2)+hj) THEN
                              i1 = i-l_start(1)-hi+1
                              j1=  j-l_start(2)-hj+1
                              OFB(i1,j1,k1,:,sample) = OFB(i1,j1,k1,:,sample) + field(i-is+1+hi,j-js+1+hj,k,:)
                        END IF
                        END DO
                     END DO
                  END DO
               ELSE IF ( reduced_k_range ) THEN
                  ksr= l_start(3)
                  ker= l_end(3)
                  OFB(is-hi:ie-hi,js-hj:je-hj,:,:,sample) = OFB(is-hi:ie-hi,js-hj:je-hj,:,:,sample) + &
                  &  field(f1:f2,f3:f4,ksr:ker,:)
               ELSE
                  IF ( debug_diag_manager ) THEN
                     CALL bbounds%update_bounds(is-hi, ie-hi, js-hj, je-hj, ks, ke)
                     CALL fms_diag_check_out_of_bounds(ofb, bbounds, output_name, module_name, err_msg_local)
                     IF ( err_msg_local /= '' ) THEN
                        IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN
                           succeded = .FALSE.
                           RETURN
                        END IF
                     END IF
                  END IF
                  OFB(is-hi:ie-hi,js-hj:je-hj,ks:ke,:,sample) = OFB(is-hi:ie-hi,js-hj:je-hj,ks:ke,:,sample) + &
                  &    field(f1:f2,f3:f4,ks:ke, :)
               END IF
            END IF MASK_PRSNT_3_IF
            count_0d = 1
            !END time_sum
         ELSE TIME_IF !! ( not average, not min, not max, not sum )
            count_0d = 1
            IF ( need_compute ) THEN
               DO j = js, je
                  DO i = is, ie
                     IF (l_start(1)+hi<= i .AND. i<= l_end(1)+hi .AND. l_start(2)+hj<= j .AND. j<= l_end(2)+hj) THEN
                        i1 = i-l_start(1)-hi+1
                        j1 = j-l_start(2)-hj+1
                        OFB(i1,j1,:,:,sample) = field(i-is+1+hi,j-js+1+hj,l_start(3):l_end(3),:)
                     END IF
                  END DO
               END DO
               ! instantaneous output
            ELSE IF ( reduced_k_range ) THEN
               ksr = l_start(3)
               ker = l_end(3)
               OFB(is-hi:ie-hi,js-hj:je-hj,:,:,sample) = field(f1:f2,f3:f4,ksr:ker,:)
            ELSE
               IF ( debug_diag_manager ) THEN
                  CALL bbounds%update_bounds(is-hi, ie-hi, js-hj, je-hj, ks, ke)
                  CALL fms_diag_check_out_of_bounds(ofb, bbounds, output_name, module_name, err_msg_local)
                  IF ( err_msg_local /= '' ) THEN
                     IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN
                        succeded = .FALSE.
                        RETURN
                     END IF
                  END IF
               END IF
               OFB(is-hi:ie-hi,js-hj:je-hj,ks:ke,:,sample) = field(f1:f2,f3:f4,ks:ke,:)
            END IF

            IF (mask_present .AND. missvalue_present ) THEN
               IF ( need_compute ) THEN
                  DO k = l_start(3), l_end(3)
                     k1 = k - l_start(3) + 1
                     DO j = js, je
                        DO i = is, ie
                           IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. &
                           & j <= l_end(2)+hj ) THEN
                              i1 = i-l_start(1)-hi+1
                              j1 =  j-l_start(2)-hj+1
                              WHERE ( .NOT.mask(i-is+1+hi,j-js+1+hj,k,:) ) &
                               & OFB(i1,j1,k1,:,sample) = missvalue
                           END IF
                        END DO
                     END DO
                  END DO
               ELSE IF ( reduced_k_range ) THEN
                  ksr= l_start(3)
                  ker= l_end(3)
                  DO k=ksr, ker
                     k1= k - ksr + 1
                     DO j=js, je
                        DO i=is, ie
                          WHERE ( mask(i-is+1+hi,j-js+1+hj,k,:) .eqv. .false.) &
                              & OFB(i-hi,j-hj,k1,:,sample)= missvalue
                        END DO
                     END DO
                  END DO
               ELSE
                  DO k=ks, ke
                     DO j=js, je
                        DO i=is, ie
                           WHERE ( .NOT. mask(i-is+1+hi,j-js+1+hj,k,:) )&
                              & OFB(i-hi,j-hj,k,:,sample)= missvalue
                        END DO
                     END DO
                  END DO
               END IF
            END IF
         END IF TIME_IF
      succeded = .TRUE.
      RETURN

   END FUNCTION fieldbuff_copy_fieldvals_r8



  !> @brief This code will be used by the  preprocessor to generate an implementation
  !! of the module procedure for the fieldbuff_copy_misvals interface. The
  !! generated function is a wrapper calling 4D field/5D buffer version of the same.
  !! TODO (MDM) the meaning of an integer rmask has to be studied.
   SUBROUTINE fieldbuff_copy_missvals_3d_r8 (ofield_cfg, ofield_index_cfg, ofb, sample, &
    & l_start, l_end, rmask, rmask_thresh, missvalue)
      TYPE(fmsDiagOutfield_type), INTENT(in) :: ofield_cfg   !< The ofield_cfg object
      TYPE(fmsDiagOutfieldIndex_type) , INTENT(in) :: ofield_index_cfg    !< The ofield_index_cfg object
      REAL(r8_kind) , CONTIGUOUS,DIMENSION(:,:,:,:),INTENT(inout),target:: ofb !<The Output Field Buffer
      INTEGER, INTENT(in) :: sample   !< index along the diurnal time axis
      REAL(r8_kind) , CONTIGUOUS,DIMENSION(:,:,:),INTENT(in),target:: rmask !<A real valued mask 4 the field
                                                      !! determined by logical mask. Updates where rmask < rmask_thresh
      REAL(r8_kind) , INTENT(in) :: rmask_thresh  !< Updates where rmask < rmask_thresh
      REAL(r8_kind) , INTENT(in) :: missvalue  !< Value used to update the buffer.
      INTEGER, DIMENSION(3), INTENT(in)  :: l_start   !< local start indices on spatial axes for regional output
      INTEGER, DIMENSION(3), INTENT(in)  :: l_end   !< local end indices on spatial for regional output

      !! These below are used in pointer bounds remapping
      REAL(r8_kind) , pointer, DIMENSION(:,:,:,:,:) :: ofb_ptr !< Pointer to the output field
                                                                        !! buffer - used in remapping
      REAL(r8_kind) , pointer, DIMENSION(:,:,:,:) :: rmask_ptr !< Pointer to the rmask - used
                                                                         !! in remapping

      !!Initialize all the pointers
      ofb_ptr(1:size(ofb,1),1:size(ofb,2),1:size(ofb,3), 1:1, 1:size(ofb,4)) => ofb
      rmask_ptr(1:size(rmask,1),1:size(rmask,2),1:size(rmask,3),1:1) => rmask

      call fieldbuff_copy_missvals_r8 (ofield_cfg, ofield_index_cfg, ofb_ptr, sample, &
        & l_start, l_end, rmask_ptr, rmask_thresh, missvalue)
  END SUBROUTINE fieldbuff_copy_missvals_3d_r8


  !> @brief This code will be used by the  preprocessor to generate an implementation
  !! of the module procedure for the fieldbuff_copy_misvals interface.
  !!  The function updates where appropriate and depending on the rmask argument,
   !! elements of the running field output buffer (argument buffer) with value missvalue.
   !! NOTE: It appears these OFB updates were introcuded by EMC MM into the tail end of the
   !! legacy send_data_3d.
   SUBROUTINE  fieldbuff_copy_missvals_r8 (ofield_cfg, ofield_index_cfg, buffer, sample, &
      & l_start, l_end, rmask, rmask_thresh, missvalue)
        TYPE(fmsDiagOutfield_type), INTENT(in) :: ofield_cfg   !< The fmsDiagOutfield_type object,
                                                               !! where "cfg" is short for configuration
        TYPE(fmsDiagOutfieldIndex_type) , INTENT(in) :: ofield_index_cfg    !< The fmsDiagOutfieldIndex_type object,
                                                               !!where "cfg" is short for configuration
        REAL(r8_kind) , INTENT(inout), DIMENSION(:,:,:,:,:) :: buffer !< the buffer to update
         INTEGER, INTENT(in) :: sample !< index along the diurnal time axis
         INTEGER, INTENT(in), DIMENSION(3):: l_start !< local start indices on 3 axes for regional output
         INTEGER, INTENT(in), DIMENSION(3):: l_end !< local end indices on 3 axes for regional output
         REAL(r8_kind) , INTENT(in), DIMENSION(:,:,:,:):: rmask  !< Updates where rmask < rmask_thresh
         REAL(r8_kind) , INTENT(in) :: rmask_thresh  !< Updates where rmask < rmask_thresh
         REAL(r8_kind) , INTENT(in) :: missvalue  !< Value used to update the buffer.

         !< Looping indices copied from corresponding one in ofield_index_cfg info:
         INTEGER :: is, js, ks, ie, je, ke, hi, hj
         !< Floags copied from corresponding one in ofield_cfg info:
         LOGICAL :: need_compute !< A copy of same variable in ofield_cfg
         LOGICAL :: reduced_k_range !< A copy of same variable in ofield_cfg
         INTEGER :: ksr, ker !< Loop indices used in reduced_k_range calculations
         !< Looping indices, derived from ofield_index_cfg info:
         INTEGER :: i, j, k, i1, j1, k1

         is = ofield_index_cfg%get_is()
         js = ofield_index_cfg%get_js()
         ks = ofield_index_cfg%get_ks()
         ie = ofield_index_cfg%get_ie()
         je = ofield_index_cfg%get_je()
         ke = ofield_index_cfg%get_ke()
         hi = ofield_index_cfg%get_hi()
         hj = ofield_index_cfg%get_hj()

         reduced_k_range = ofield_cfg%get_reduced_k_range()
         need_compute = ofield_cfg%get_need_compute()

         associate(ofb => buffer)

            ! If rmask and missing value present, then insert missing value
            IF ( need_compute ) THEN
               DO k = l_start(3), l_end(3)
                  k1 = k - l_start(3) + 1
                  DO j = js, je
                     DO i = is, ie
                        IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND.&
                        & j <= l_end(2)+hj ) THEN
                           i1 = i-l_start(1)-hi+1
                           j1 =  j-l_start(2)-hj+1
                           where ( rmask(i-is+1+hi,j-js+1+hj,k,:) <= rmask_thresh )
                              ofb(i1,j1,k1,:,sample) = missvalue
                           end where
                        END IF
                     END DO
                  END DO
               END DO
            ELSE IF ( reduced_k_range ) THEN
               ksr= l_start(3)
               ker= l_end(3)
               DO k= ksr, ker
                  k1 = k - ksr + 1
                  DO j=js, je
                     DO i=is, ie
                        where ( rmask(i-is+1+hi,j-js+1+hj,k,:) <= rmask_thresh )
                           ofb(i-hi,j-hj,k1,:,sample)= missvalue
                        endwhere
                     END DO
                  END DO
               END DO
            ELSE
               DO k=ks, ke
                  DO j=js, je
                     DO i=is, ie
                        where ( rmask(i-is+1+hi,j-js+1+hj,k,:) <= rmask_thresh )
                          ofb(i-hi,j-hj,k,:,sample)= missvalue
                        endwhere
                     END DO
                  END DO
               END DO
            END IF
         end associate
      END SUBROUTINE  fieldbuff_copy_missvals_r8
# 50 "diag_manager/include/fms_diag_fieldbuff_update.inc" 2
# 98 "diag_manager/fms_diag_fieldbuff_update.F90" 2

END MODULE fms_diag_fieldbuff_update_mod
!> @}
! close documentation grouping
