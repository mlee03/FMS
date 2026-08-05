# 1 "monin_obukhov/monin_obukhov.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "monin_obukhov/monin_obukhov.F90"
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
!> @defgroup monin_obukhov_mod monin_obukhov_mod
!! @ingroup monin_obukhov
!! @{
!! @brief Routines for computing surface drag coefficients
!! from data at the lowest model level
!! and for computing the profile of fields
!! between the lowest model level and the ground
!! using Monin-Obukhov scaling

module monin_obukhov_mod

use constants_mod, only: grav, vonkarm
use mpp_mod,       only: input_nml_file
use fms_mod,       only: error_mesg, FATAL, check_nml_error,   &
                         mpp_pe, mpp_root_pe, stdlog, &
                         write_version_number
use monin_obukhov_inter, only: monin_obukhov_diff, monin_obukhov_drag_1d, &
                               monin_obukhov_profile_1d, monin_obukhov_stable_mix
use platform_mod,        only: r4_kind, r8_kind
implicit none
private

!=======================================================================
 public :: monin_obukhov_init
 public :: monin_obukhov_end
 public :: mo_drag
 public :: mo_profile
 public :: mo_diff
 public :: stable_mix
!=======================================================================

!> @brief Compute surface drag coefficients
interface mo_drag
    module procedure mo_drag_0d_r4, mo_drag_0d_r8
    module procedure mo_drag_1d_r4, mo_drag_1d_r8
    module procedure mo_drag_2d_r4, mo_drag_2d_r8
end interface


interface mo_profile
    module procedure mo_profile_0d_r4, mo_profile_0d_r8
    module procedure mo_profile_1d_r4, mo_profile_1d_r8
    module procedure mo_profile_2d_r4, mo_profile_2d_r8
    module procedure mo_profile_0d_n_r4, mo_profile_0d_n_r8
    module procedure mo_profile_1d_n_r4, mo_profile_1d_n_r8
    module procedure mo_profile_2d_n_r4, mo_profile_2d_n_r8
end interface

interface mo_diff
    module procedure mo_diff_0d_n_r4, mo_diff_0d_n_r8
    module procedure mo_diff_0d_1_r4, mo_diff_0d_1_r8
    module procedure mo_diff_1d_n_r4, mo_diff_1d_n_r8
    module procedure mo_diff_1d_1_r4, mo_diff_1d_1_r8
    module procedure mo_diff_2d_n_r4, mo_diff_2d_n_r8
    module procedure mo_diff_2d_1_r4, mo_diff_2d_1_r8
end interface

interface stable_mix
    module procedure stable_mix_0d_r4, stable_mix_0d_r8
    module procedure stable_mix_1d_r4, stable_mix_1d_r8
    module procedure stable_mix_2d_r4, stable_mix_2d_r8
    module procedure stable_mix_3d_r4, stable_mix_3d_r8
end interface

interface mo_integral_m
    module procedure mo_integral_m_r4, mo_integral_m_r8
end interface mo_integral_m

interface mo_integral_tq
    module procedure mo_integral_tq_r4, mo_integral_tq_r8
end interface mo_integral_tq

interface mo_derivative_m
    module procedure mo_derivative_m_r4, mo_derivative_m_r8
end interface mo_derivative_m

interface mo_derivative_t
    module procedure mo_derivative_t_r4, mo_derivative_t_r8
end interface mo_derivative_t

!-----------------------------------------------------------------------
! version number of this module
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
# 102 "monin_obukhov/monin_obukhov.F90" 2

!=======================================================================

!  DEFAULT VALUES OF NAMELIST PARAMETERS:

real(kind=r8_kind) :: rich_crit      = 2.0_r8_kind
real(kind=r8_kind) :: drag_min_heat  = 1.0E-05_r8_kind
real(kind=r8_kind) :: drag_min_moist = 1.0E-05_r8_kind
real(kind=r8_kind) :: drag_min_mom   = 1.0E-05_r8_kind
logical            :: neutral        = .false.
integer            :: stable_option  = 1
real(kind=r8_kind) :: zeta_trans     = 0.5_r8_kind
logical            :: new_mo_option  = .false.


namelist /monin_obukhov_nml/ rich_crit, neutral, drag_min_heat, &
                             drag_min_moist, drag_min_mom,      &
                             stable_option, zeta_trans, new_mo_option !miz

!=======================================================================

!  MODULE VARIABLES

real(kind=r8_kind), parameter    :: small  = 1.0E-04_r8_kind
real(kind=r8_kind)               :: b_stab, r_crit, lambda, rich_trans
real(kind=r8_kind)               :: sqrt_drag_min_heat, sqrt_drag_min_moist, sqrt_drag_min_mom
logical                          :: module_is_initialized = .false.


contains

!=======================================================================

subroutine monin_obukhov_init

integer :: ierr, io, logunit

!------------------- read namelist input -------------------------------

      read (input_nml_file, nml=monin_obukhov_nml, iostat=io)
      ierr = check_nml_error(io,"monin_obukhov_nml")

!---------- output namelist to log-------------------------------------

      if ( mpp_pe() == mpp_root_pe() ) then
           call write_version_number('MONIN_OBUKOV_MOD', version)
           logunit = stdlog()
           write (logunit, nml=monin_obukhov_nml)
      endif

!----------------------------------------------------------------------

if(rich_crit.le.0.25_r8_kind)  call error_mesg( &
        'MONIN_OBUKHOV_INIT in MONIN_OBUKHOV_MOD', &
        'rich_crit in monin_obukhov_mod must be > 0.25', FATAL)

if(drag_min_heat.le.0.0_r8_kind)  call error_mesg( &
        'MONIN_OBUKHOV_INIT in MONIN_OBUKHOV_MOD', &
        'drag_min_heat in monin_obukhov_mod must be >= 0.0', FATAL)

if(drag_min_moist.le.0.0_r8_kind)  call error_mesg( &
        'MONIN_OBUKHOV_INIT in MONIN_OBUKHOV_MOD', &
        'drag_min_moist in monin_obukhov_mod must be >= 0.0', FATAL)

if(drag_min_mom.le.0.0_r8_kind)  call error_mesg( &
        'MONIN_OBUKHOV_INIT in MONIN_OBUKHOV_MOD', &
        'drag_min_mom in monin_obukhov_mod must be >= 0.0', FATAL)

if(stable_option < 1 .or. stable_option > 2) call error_mesg( &
        'MONIN_OBUKHOV_INIT in MONIN_OBUKHOV_MOD', &
        'the only allowable values of stable_option are 1 and 2', FATAL)

if(stable_option == 2 .and. zeta_trans < 0) call error_mesg( &
        'MONIN_OBUKHOV_INIT in MONIN_OBUKHOV_MOD', &
        'zeta_trans must be positive', FATAL)

b_stab = 1.0_r8_kind/rich_crit
r_crit = 0.95_r8_kind*rich_crit  ! convergence can get slow if one is
                         ! close to rich_crit

sqrt_drag_min_heat = 0.0_r8_kind
if(drag_min_heat.ne.0.0_r8_kind) sqrt_drag_min_heat = sqrt(drag_min_heat)

sqrt_drag_min_moist = 0.0_r8_kind
if(drag_min_moist.ne.0.0_r8_kind) sqrt_drag_min_moist = sqrt(drag_min_moist)

sqrt_drag_min_mom = 0.0_r8_kind
if(drag_min_mom.ne.0.0_r8_kind) sqrt_drag_min_mom = sqrt(drag_min_mom)

lambda     = 1.0_r8_kind + (5.0_r8_kind - b_stab)*zeta_trans   ! used only if stable_option = 2
rich_trans = zeta_trans/(1.0_r8_kind + 5.0_r8_kind*zeta_trans) ! used only if stable_option = 2

module_is_initialized = .true.

return
end subroutine monin_obukhov_init

!=======================================================================

subroutine monin_obukhov_end

module_is_initialized = .false.

end subroutine monin_obukhov_end

!=======================================================================


# 1 "monin_obukhov/include/monin_obukhov_r4.fh" 1
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












































































# 1 "monin_obukhov/include/monin_obukhov.inc" 1
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


subroutine mo_drag_1d_r4 &
         (pt, pt0, z, z0, zt, zq, speed, drag_m, drag_t, drag_q, &
          u_star, b_star, avail)

real(kind=r4_kind), intent(in)   , dimension(:) :: pt, pt0, z, z0, zt, zq, speed
real(kind=r4_kind), intent(inout), dimension(:) :: drag_m, drag_t, drag_q, u_star, b_star
logical, intent(in), optional, dimension(:)          :: avail

logical                                              :: lavail, avail_dummy(1)
integer                                              :: n, ier

integer, parameter                                   :: max_iter = 20
integer, parameter                                   :: lkind    = r4_kind
real(kind=r4_kind), parameter                   :: error    = 1.0E-04_lkind, &
                                                        zeta_min = 1.0E-06_lkind, &
                                                        small    = 1.0E-04_lkind

! #include "monin_obukhov_interfaces.h"

if(.not.module_is_initialized) call error_mesg('mo_drag_1d in monin_obukhov_mod', &
     'monin_obukhov_init has not been called', FATAL)

n      = size(pt)
lavail = .false.
if(present(avail)) lavail = .true.


if(lavail) then
   if (count(avail) .eq. 0) return
   call monin_obukhov_drag_1d(real(grav, r4_kind), real(vonkarm, r4_kind),    &
        & error, zeta_min, max_iter, real(small, r4_kind), neutral, stable_option, &
        & new_mo_option, real(rich_crit, r4_kind), real(zeta_trans, r4_kind), &!miz
        & real(drag_min_heat, r4_kind), real(drag_min_moist, r4_kind),        &
        & real(drag_min_mom, r4_kind), n, pt, pt0, z, z0, zt,                      &
        & zq, speed, drag_m, drag_t, drag_q, u_star, b_star, lavail, avail, ier)
else
call monin_obukhov_drag_1d(real(grav, r4_kind), real(vonkarm, r4_kind),             &
        & error, zeta_min, max_iter, real(small, r4_kind), neutral, stable_option,       &
        & new_mo_option, real(rich_crit, r4_kind), real(zeta_trans, r4_kind),       &!miz
        & real(drag_min_heat, r4_kind), real(drag_min_moist, r4_kind),              &
        & real(drag_min_mom, r4_kind), n, pt, pt0, z, z0, zt, zq, speed, drag_m, drag_t, &
        & drag_q, u_star, b_star, lavail, avail_dummy, ier)
endif

end subroutine mo_drag_1d_r4


!=======================================================================

subroutine mo_profile_1d_r4(zref, zref_t, z, z0, zt, zq, u_star, b_star, q_star, &
                         del_m, del_t, del_q, avail)

real(kind=r4_kind),    intent(in)                :: zref, zref_t
real(kind=r4_kind),    intent(in) , dimension(:) :: z, z0, zt, zq, u_star, b_star, q_star
real(kind=r4_kind),    intent(out), dimension(:) :: del_m, del_t, del_q
logical, intent(in) , optional, dimension(:)          :: avail

logical                                               :: dummy_avail(1)
integer                                               :: n, ier

! #include "monin_obukhov_interfaces.h"

if(.not.module_is_initialized) call error_mesg('mo_profile_1d in monin_obukhov_mod', &
     'monin_obukhov_init has not been called', FATAL)

n = size(z)
if(present(avail)) then

   if (count(avail) .eq. 0) return

   call monin_obukhov_profile_1d(real(vonkarm, r4_kind), &
        & neutral, stable_option, new_mo_option, real(rich_crit, r4_kind),   &
        & real(zeta_trans, r4_kind), n, zref, zref_t, z, z0, zt, zq, u_star, &
        & b_star, q_star, del_m, del_t, del_q, .true., avail, ier)

else

   call monin_obukhov_profile_1d(real(vonkarm, r4_kind), &
        & neutral, stable_option, new_mo_option, real(rich_crit, r4_kind),   &
        & real(zeta_trans, r4_kind), n, zref, zref_t, z, z0, zt, zq, u_star, &
        & b_star, q_star, del_m, del_t, del_q, .false., dummy_avail, ier)

endif

end subroutine mo_profile_1d_r4

!=======================================================================

subroutine stable_mix_3d_r4(rich, mix)

real(kind=r4_kind), intent(in) , dimension(:,:,:)  :: rich
real(kind=r4_kind), intent(out), dimension(:,:,:)  :: mix
integer :: n2 !< Size of dimension 2 of mix and rich
integer :: n3 !< Size of dimension 3 of mix and rich
integer :: i, j !< Loop indices

n2 = size(mix, 2)
n3 = size(mix, 3)

do j=1, n3
  do i=1, n2
    call stable_mix(rich(:, i, j), mix(:, i, j))
  enddo
enddo

end subroutine stable_mix_3d_r4

!=======================================================================

subroutine mo_diff_2d_n_r4(z, u_star, b_star, k_m, k_h)

real(kind=r4_kind), intent(in),  dimension(:,:,:) :: z
real(kind=r4_kind), intent(in),  dimension(:,:)   :: u_star, b_star
real(kind=r4_kind), intent(out), dimension(:,:,:) :: k_m, k_h

integer                                                :: ni, nj, nk, ier
integer, parameter                                     :: lkind = r4_kind
real(kind=r4_kind), parameter                     :: ustar_min = 1.0E-10_lkind

if(.not.module_is_initialized) call error_mesg('mo_diff_2d_n in monin_obukhov_mod', &
     'monin_obukhov_init has not been called', FATAL)

ni = size(z, 1); nj = size(z, 2); nk = size(z, 3)
call monin_obukhov_diff(real(vonkarm, r4_kind), ustar_min, neutral,                      &
                             & stable_option, new_mo_option, real(rich_crit, r4_kind),   &
                             & real(zeta_trans, r4_kind), ni, nj, nk, z, u_star, b_star, &
                             & k_m, k_h, ier)

end subroutine mo_diff_2d_n_r4

!=======================================================================
! The following routines are used by the public interfaces above
!=======================================================================

subroutine solve_zeta_r4(rich, z, z0, zt, zq, f_m, f_t, f_q, mask)

real(kind=r4_kind), intent(in) , dimension(:) :: rich, z, z0, zt, zq
logical, intent(in) , dimension(:)                 :: mask
real(kind=r4_kind), intent(out), dimension(:) :: f_m, f_t, f_q

integer, parameter                                 :: lkind    = r4_kind
real(kind=r4_kind), parameter                 :: error    = 1.0E-04_lkind
real(kind=r4_kind), parameter                 :: zeta_min = 1.0E-06_lkind
integer, parameter                                 :: max_iter = 20

real(kind=r4_kind)                            :: max_cor
integer                                            :: iter

real(kind=r4_kind), dimension(size(rich(:)))  ::   &
          d_rich, rich_1, correction, corr, z_z0, z_zt, z_zq, &
          ln_z_z0, ln_z_zt, ln_z_zq, zeta,                    &
          phi_m, phi_m_0, phi_t, phi_t_0, rzeta,              &
          zeta_0, zeta_t, zeta_q, df_m, df_t

logical, dimension(size(rich(:)))                  :: mask_1

z_z0    = z/z0
z_zt    = z/zt
z_zq    = z/zq
ln_z_z0 = log(z_z0)
ln_z_zt = log(z_zt)
ln_z_zq = log(z_zq)

corr = 0.0_lkind
mask_1 = mask

! initial guess

where(mask_1)
  zeta = rich*ln_z_z0*ln_z_z0/ln_z_zt
elsewhere
  zeta = 0.0_lkind
end where

where (mask_1 .and. rich >= 0.0_lkind)
  zeta = zeta/(1.0_lkind - rich/real(rich_crit, r4_kind))
end where

iter_loop: do iter = 1, max_iter

  where (mask_1 .and. abs(zeta).lt.zeta_min)
    zeta   = 0.0_lkind
    f_m    = ln_z_z0
    f_t    = ln_z_zt
    f_q    = ln_z_zq
    mask_1 = .false.  ! don't do any more calculations at these pts
  end where

  where (mask_1)
    rzeta  = 1.0_lkind/zeta
    zeta_0 = zeta/z_z0
    zeta_t = zeta/z_zt
    zeta_q = zeta/z_zq
  elsewhere
    zeta_0 = 0.0_lkind
    zeta_t = 0.0_lkind
    zeta_q = 0.0_lkind
  end where

  call mo_derivative_m(phi_m  , zeta  , mask_1)
  call mo_derivative_m(phi_m_0, zeta_0, mask_1)
  call mo_derivative_t(phi_t  , zeta  , mask_1)
  call mo_derivative_t(phi_t_0, zeta_t, mask_1)

  call mo_integral_m(f_m, zeta, zeta_0, ln_z_z0, mask_1)
  call mo_integral_tq(f_t, f_q, zeta, zeta_t, zeta_q, ln_z_zt, ln_z_zq, mask_1)

  where (mask_1)
    df_m       = (phi_m - phi_m_0)*rzeta
    df_t       = (phi_t - phi_t_0)*rzeta
    rich_1     = zeta*f_t/(f_m*f_m)
    d_rich     = rich_1*( rzeta +  df_t/f_t - 2.0_lkind *df_m/f_m)
    correction = (rich - rich_1)/d_rich
    corr       = min(abs(correction),abs(correction/zeta))
      ! the criterion corr < error seems to work ok, but is a bit arbitrary
      !  when zeta is small the tolerance is reduced
  end where

  max_cor = maxval(corr)

  if(max_cor > error) then
    mask_1 = mask_1 .and. (corr > error)
       ! change the mask so computation proceeds only on non-converged points
    where(mask_1)
      zeta = zeta + correction
    end where
    cycle iter_loop
  else
    return
  end if

end do iter_loop

call error_mesg ('solve_zeta in monin_obukhov_mod',  &
                 'surface drag iteration did not converge', FATAL)

end subroutine solve_zeta_r4

!=======================================================================

subroutine mo_derivative_m_r4(phi_m, zeta, mask)

! the differential similarity function for momentum

real(kind=r4_kind), intent(out),  dimension(:) :: phi_m
real(kind=r4_kind), intent(in),   dimension(:) :: zeta
logical,                 intent(in),   dimension(:) :: mask

logical, dimension(size(zeta(:)))                   :: stable, unstable
real(kind=r4_kind), dimension(size(zeta(:)))   :: x
integer, parameter                                  :: lkind = r4_kind

stable   = mask .and. zeta >= 0.0_lkind
unstable = mask .and. zeta <  0.0_lkind

where (unstable)
  x     = (1.0_lkind - 16.0_lkind*zeta  )**(-0.5_lkind)
  phi_m = sqrt(x)  ! phi_m = (1 - 16.0*zeta)**(-0.25)
end where

if(stable_option == 1) then

  where (stable)
    phi_m = 1.0_lkind + zeta*(5.0_lkind + real(b_stab, r4_kind) &
            *zeta)/(1.0_lkind + zeta)
  end where

else if(stable_option == 2) then

  where (stable .and. zeta < real(zeta_trans,r4_kind))
    phi_m = 1.0_lkind + 5.0_lkind*zeta
  end where
  where (stable .and. zeta >= real(zeta_trans,r4_kind))
    phi_m = real(lambda, r4_kind) + real(b_stab, r4_kind)*zeta
  end where

endif

return
end subroutine mo_derivative_m_r4

!=======================================================================

subroutine mo_derivative_t_r4(phi_t, zeta, mask)

! the differential similarity function for buoyancy and tracers

real(kind=r4_kind), intent(out),  dimension(:) :: phi_t
real(kind=r4_kind), intent(in),   dimension(:) :: zeta
logical ,                intent(in),   dimension(:) :: mask

logical, dimension(size(zeta(:)))                   :: stable, unstable
integer, parameter                                  :: lkind = r4_kind

stable   = mask .and. zeta >= 0.0_lkind
unstable = mask .and. zeta <  0.0_lkind

where (unstable)
  phi_t = (1.0_lkind - 16.0_lkind*zeta)**(-0.5_lkind)
end where

if(stable_option == 1) then

  where (stable)
    phi_t = 1.0_lkind + zeta * (5.0_lkind + real(b_stab, r4_kind)&
            * zeta)/(1.0_lkind + zeta)
  end where

else if(stable_option == 2) then

  where (stable .and. zeta < real(zeta_trans,r4_kind))
    phi_t = 1.0_lkind + 5.0_lkind*zeta
  end where
  where (stable .and. zeta >= real(zeta_trans,r4_kind))
    phi_t = real(lambda, r4_kind) + real(b_stab, r4_kind)*zeta
  end where

endif

return
end subroutine mo_derivative_t_r4

!=======================================================================

subroutine mo_integral_tq_r4 (psi_t, psi_q, zeta, zeta_t, zeta_q, &
                           ln_z_zt, ln_z_zq, mask)

! the integral similarity function for moisture and tracers

real(kind=r4_kind), intent(out), dimension(:) :: psi_t, psi_q
real(kind=r4_kind), intent(in),  dimension(:) :: zeta, zeta_t, zeta_q, ln_z_zt, ln_z_zq
logical , intent(in),  dimension(:)                :: mask

real(kind=r4_kind), dimension(size(zeta(:)))  :: x, x_t, x_q

logical, dimension(size(zeta(:)))                  :: stable, unstable, &
                                                      weakly_stable, strongly_stable
integer, parameter                                 :: lkind = r4_kind

stable   = mask .and. zeta >= 0.0_lkind
unstable = mask .and. zeta <  0.0_lkind

where(unstable)

  x     = sqrt(1.0_lkind - 16.0_lkind*zeta)
  x_t   = sqrt(1.0_lkind - 16.0_lkind*zeta_t)
  x_q   = sqrt(1.0_lkind - 16.0_lkind*zeta_q)

  psi_t = ln_z_zt - 2.0_lkind*log( (1.0_lkind + x)/(1.0_lkind + x_t) )
  psi_q = ln_z_zq - 2.0_lkind*log( (1.0_lkind + x)/(1.0_lkind + x_q) )

end where

if( stable_option == 1) then

  where (stable)

    psi_t = ln_z_zt + (5.0_lkind - real(b_stab, r4_kind)) &
            *log((1.0_lkind + zeta)/(1.0_lkind + zeta_t)) &
            + real(b_stab, r4_kind)*(zeta - zeta_t)
    psi_q = ln_z_zq + (5.0_lkind - real(b_stab, r4_kind)) &
            *log((1.0_lkind + zeta)/(1.0_lkind + zeta_q)) &
            + real(b_stab, r4_kind)*(zeta - zeta_q)

  end where

else if (stable_option == 2) then

  weakly_stable   = stable .and. zeta <= real(zeta_trans,r4_kind)
  strongly_stable = stable .and. zeta >  real(zeta_trans,r4_kind)

  where (weakly_stable)
    psi_t = ln_z_zt + 5.0_lkind*(zeta - zeta_t)
    psi_q = ln_z_zq + 5.0_lkind*(zeta - zeta_q)
  end where

  where(strongly_stable)
    x = (real(lambda, r4_kind) - 1.0_lkind)*log(zeta/real(zeta_trans, r4_kind)) + &
         real(b_stab, r4_kind)*(zeta - real(zeta_trans, r4_kind))
  endwhere

  where (strongly_stable .and. zeta_t <= real(zeta_trans,r4_kind))
    psi_t = ln_z_zt + x + 5.0_lkind * (real(zeta_trans, r4_kind) - zeta_t)
  end where

  where (strongly_stable .and. zeta_t > real(zeta_trans,r4_kind))
    psi_t = real(lambda, r4_kind)* ln_z_zt &
            + real(b_stab, r4_kind)*(zeta - zeta_t)
  endwhere

  where (strongly_stable .and. zeta_q <= real(zeta_trans,r4_kind))
    psi_q = ln_z_zq + x + 5.0_lkind &
            *(real(zeta_trans, r4_kind) - zeta_q)
  end where

  where (strongly_stable .and. zeta_q > real(zeta_trans,r4_kind))
    psi_q = real(lambda, r4_kind)*ln_z_zq + real(b_stab, r4_kind) &
            * (zeta - zeta_q)
  endwhere

end if

return
end subroutine mo_integral_tq_r4

!=======================================================================

subroutine mo_integral_m_r4 (psi_m, zeta, zeta_0, ln_z_z0, mask)

!  the integral similarity function for momentum

real(kind=r4_kind), intent(out), dimension(:) :: psi_m
real(kind=r4_kind), intent(in),  dimension(:) :: zeta, zeta_0, ln_z_z0
logical,                 intent(in),  dimension(:) :: mask

real(kind=r4_kind), dimension(size(zeta(:)))  :: x, x_0, x1, x1_0, num, denom, y

logical, dimension(size(zeta(:)))                  :: stable, unstable, &
                                                      weakly_stable, strongly_stable
integer, parameter                                 :: lkind = r4_kind

stable   = mask .and. zeta >= 0.0_lkind
unstable = mask .and. zeta <  0.0_lkind

where(unstable)

  x      = sqrt(1.0_lkind - 16.0_lkind*zeta)
  x_0    = sqrt(1.0_lkind - 16.0_lkind*zeta_0)

  x      = sqrt(x)
  x_0    = sqrt(x_0)

  x1     = 1.0_lkind + x
  x1_0   = 1.0_lkind + x_0

  num    = x1*x1*(1.0_lkind + x*x)
  denom  = x1_0*x1_0*(1.0_lkind + x_0*x_0)
  y      = atan(x) - atan(x_0)
  psi_m  = ln_z_z0 - log(num/denom) + 2.0_lkind*y

end where

if( stable_option == 1) then

  where (stable)
    psi_m = ln_z_z0 + (5.0_lkind - real(b_stab, r4_kind)) &
            *log((1.0_lkind + zeta)/(1.0_lkind + zeta_0)) &
            + real(b_stab, r4_kind)*(zeta - zeta_0)
  end where

else if (stable_option == 2) then

  weakly_stable   = stable .and. zeta <= real(zeta_trans,r4_kind)
  strongly_stable = stable .and. zeta >  real(zeta_trans,r4_kind)

  where (weakly_stable)
    psi_m = ln_z_z0 + 5.0_lkind*(zeta - zeta_0)
  end where

  where(strongly_stable)
    x = (real(lambda, r4_kind) - 1.0_lkind)*log(zeta/real(zeta_trans, r4_kind)) + &
        real(b_stab, r4_kind)*(zeta - real(zeta_trans, r4_kind))
  endwhere

  where (strongly_stable .and. zeta_0 <= real(zeta_trans,r4_kind))
    psi_m = ln_z_z0 + x + 5.0_lkind &
            *(real(zeta_trans, r4_kind) - zeta_0)
  end where
  where (strongly_stable .and. zeta_0 > real(zeta_trans,r4_kind))
    psi_m = real(lambda, r4_kind)*ln_z_z0 + real(b_stab, r4_kind) &
            *(zeta - zeta_0)
  endwhere

end if

return
end subroutine mo_integral_m_r4


!=======================================================================
! The following routines allow the public interfaces to be used
! with different dimensions of the input and output
!
!=======================================================================


subroutine mo_drag_2d_r4 &
    (pt, pt0, z, z0, zt, zq, speed, drag_m, drag_t, drag_q, u_star, b_star)

real(kind=r4_kind), intent(in)   , dimension(:,:) :: z, speed, pt, pt0, z0, zt, zq
real(kind=r4_kind), intent(out)  , dimension(:,:) :: drag_m, drag_t, drag_q
real(kind=r4_kind), intent(inout), dimension(:,:) :: u_star, b_star

integer :: j

do j = 1, size(pt,2)
  call mo_drag (pt(:,j), pt0(:,j), z(:,j), z0(:,j), zt(:,j), zq(:,j), &
                   speed(:,j), drag_m(:,j), drag_t(:,j), drag_q(:,j), &
                   u_star(:,j), b_star(:,j))
end do


return
end subroutine mo_drag_2d_r4

!=======================================================================
subroutine mo_drag_0d_r4 &
    (pt, pt0, z, z0, zt, zq, speed, drag_m, drag_t, drag_q, u_star, b_star)

real(kind=r4_kind), intent(in)    :: z, speed, pt, pt0, z0, zt, zq
real(kind=r4_kind), intent(out)   :: drag_m, drag_t, drag_q, u_star, b_star

real(kind=r4_kind), dimension(1)  :: pt_1, pt0_1, z_1, z0_1, zt_1, zq_1, speed_1, &
                      drag_m_1, drag_t_1, drag_q_1, u_star_1, b_star_1

pt_1   (1) = pt
pt0_1  (1) = pt0
z_1    (1) = z
z0_1   (1) = z0
zt_1   (1) = zt
zq_1   (1) = zq
speed_1(1) = speed

call mo_drag (pt_1, pt0_1, z_1, z0_1, zt_1, zq_1, speed_1, &
                 drag_m_1, drag_t_1, drag_q_1, u_star_1, b_star_1)

drag_m = drag_m_1(1)
drag_t = drag_t_1(1)
drag_q = drag_q_1(1)
u_star = u_star_1(1)
b_star = b_star_1(1)

return
end subroutine mo_drag_0d_r4
!=======================================================================

subroutine mo_profile_2d_r4(zref, zref_t, z, z0, zt, zq, u_star, b_star, q_star, &
                         del_m, del_h, del_q)

real(kind=r4_kind), intent(in)                  :: zref, zref_t
real(kind=r4_kind), intent(in) , dimension(:,:) :: z, z0, zt, zq, u_star, b_star, q_star
real(kind=r4_kind), intent(out), dimension(:,:) :: del_m, del_h, del_q

integer :: j

do j = 1, size(z,2)
  call mo_profile (zref, zref_t, z(:,j), z0(:,j), zt(:,j),         &
                      zq(:,j), u_star(:,j), b_star(:,j), q_star(:,j), &
                      del_m(:,j), del_h (:,j), del_q (:,j))
enddo

return
end subroutine mo_profile_2d_r4

!=======================================================================

subroutine mo_profile_0d_r4(zref, zref_t, z, z0, zt, zq, u_star, b_star, q_star, &
                         del_m, del_h, del_q)

real(kind=r4_kind), intent(in)  :: zref, zref_t
real(kind=r4_kind), intent(in)  :: z, z0, zt, zq, u_star, b_star, q_star
real(kind=r4_kind), intent(out) :: del_m, del_h, del_q

real(kind=r4_kind), dimension(1) :: z_1, z0_1, zt_1, zq_1, u_star_1, b_star_1, q_star_1, &
                      del_m_1, del_h_1, del_q_1

z_1     (1) = z
z0_1    (1) = z0
zt_1    (1) = zt
zq_1    (1) = zq
u_star_1(1) = u_star
b_star_1(1) = b_star
q_star_1(1) = q_star

call mo_profile (zref, zref_t, z_1, z0_1, zt_1, zq_1, &
                    u_star_1, b_star_1, q_star_1,        &
                    del_m_1, del_h_1, del_q_1)

del_m = del_m_1(1)
del_h = del_h_1(1)
del_q = del_q_1(1)


return
end subroutine mo_profile_0d_r4

!=======================================================================

subroutine mo_profile_1d_n_r4(zref, z, z0, zt, zq, u_star, b_star, q_star, &
                         del_m, del_t, del_q, avail)

real(kind=r4_kind),    intent(in),  dimension(:)   :: zref
real(kind=r4_kind),    intent(in) , dimension(:)   :: z, z0, zt, zq, u_star, b_star, q_star
real(kind=r4_kind),    intent(out), dimension(:,:) :: del_m, del_t, del_q
logical, intent(in) , optional,          dimension(:)   :: avail

integer :: k

do k = 1, size(zref(:))
  if(present(avail)) then
    call mo_profile (zref(k), zref(k), z, z0, zt, zq, &
       u_star, b_star, q_star, del_m(:,k), del_t(:,k), del_q(:,k), avail)
  else
      call mo_profile (zref(k), zref(k), z, z0, zt, zq, &
       u_star, b_star, q_star, del_m(:,k), del_t(:,k), del_q(:,k))
  endif
enddo

return
end subroutine mo_profile_1d_n_r4

!=======================================================================

subroutine mo_profile_0d_n_r4(zref, z, z0, zt, zq, u_star, b_star, q_star, &
                         del_m, del_t, del_q)

real(kind=r4_kind),    intent(in),  dimension(:) :: zref
real(kind=r4_kind),    intent(in)                :: z, z0, zt, zq, u_star, b_star, q_star
real(kind=r4_kind),    intent(out), dimension(:) :: del_m, del_t, del_q

integer :: k

do k = 1, size(zref(:))
  call mo_profile (zref(k), zref(k), z, z0, zt, zq, &
       u_star, b_star, q_star, del_m(k), del_t(k), del_q(k))
enddo

return
end subroutine mo_profile_0d_n_r4

!=======================================================================

subroutine mo_profile_2d_n_r4(zref, z, z0, zt, zq, u_star, b_star, q_star, &
                         del_m, del_t, del_q)

real(kind=r4_kind),    intent(in),  dimension(:)     :: zref
real(kind=r4_kind),    intent(in),  dimension(:,:)   :: z, z0, zt, zq, u_star, b_star, q_star
real(kind=r4_kind),    intent(out), dimension(:,:,:) :: del_m, del_t, del_q

integer :: k

do k = 1, size(zref(:))
  call mo_profile (zref(k), zref(k), z, z0, zt, zq, &
       u_star, b_star, q_star, del_m(:,:,k), del_t(:,:,k), del_q(:,:,k))
enddo

return
end subroutine mo_profile_2d_n_r4

!=======================================================================

subroutine mo_diff_2d_1_r4(z, u_star, b_star, k_m, k_h)

real(kind=r4_kind), intent(in),  dimension(:,:)      :: z, u_star, b_star
real(kind=r4_kind), intent(out), dimension(:,:)      :: k_m, k_h

real(kind=r4_kind), dimension(size(z,1),size(z,2),1) :: z_n, k_m_n, k_h_n

z_n(:,:,1) = z

call mo_diff(z_n, u_star, b_star, k_m_n, k_h_n)

k_m = k_m_n(:,:,1)
k_h = k_h_n(:,:,1)

return
end subroutine mo_diff_2d_1_r4


!=======================================================================

subroutine mo_diff_1d_1_r4(z, u_star, b_star, k_m, k_h)

real(kind=r4_kind), intent(in),  dimension(:) :: z, u_star, b_star
real(kind=r4_kind), intent(out), dimension(:) :: k_m, k_h

real(kind=r4_kind), dimension(size(z),1,1)    :: z_n, k_m_n, k_h_n
real(kind=r4_kind), dimension(size(z),1)      :: u_star_n, b_star_n

z_n   (:,1,1) = z
u_star_n(:,1) = u_star
b_star_n(:,1) = b_star

call mo_diff(z_n, u_star_n, b_star_n, k_m_n, k_h_n)

k_m = k_m_n(:,1,1)
k_h = k_h_n(:,1,1)

return
end subroutine mo_diff_1d_1_r4

!=======================================================================

subroutine mo_diff_1d_n_r4(z, u_star, b_star, k_m, k_h)

real(kind=r4_kind), intent(in),  dimension(:,:) :: z
real(kind=r4_kind), intent(in),  dimension(:)   :: u_star, b_star
real(kind=r4_kind), intent(out), dimension(:,:) :: k_m, k_h

real(kind=r4_kind), dimension(size(z,1),1)            :: u_star2, b_star2
real(kind=r4_kind), dimension(size(z,1),1, size(z,2)) :: z2, k_m2, k_h2

integer :: n

do n = 1, size(z,2)
  z2   (:,1,n) = z(:,n)
enddo
u_star2(:,1) = u_star
b_star2(:,1) = b_star

call mo_diff(z2, u_star2, b_star2, k_m2, k_h2)

do n = 1, size(z,2)
  k_m(:,n) = k_m2(:,1,n)
  k_h(:,n) = k_h2(:,1,n)
enddo

return
end subroutine mo_diff_1d_n_r4

!=======================================================================

subroutine mo_diff_0d_1_r4(z, u_star, b_star, k_m, k_h)

real(kind=r4_kind), intent(in)       :: z, u_star, b_star
real(kind=r4_kind), intent(out)      :: k_m, k_h

integer                                   :: ni, nj, nk, ier
integer, parameter                        :: lkind = r4_kind
real(kind=r4_kind), parameter        :: ustar_min = 1.0E-10_lkind
real(kind=r4_kind), dimension(1,1,1) :: z_a, k_m_a, k_h_a
real(kind=r4_kind), dimension(1,1)   :: u_star_a, b_star_a

if(.not.module_is_initialized) call error_mesg('mo_diff_0d_1 in monin_obukhov_mod', &
     'monin_obukhov_init has not been called', FATAL)

ni = 1; nj = 1; nk = 1
z_a(1,1,1)    = z
u_star_a(1,1) = u_star
b_star_a(1,1) = b_star
call monin_obukhov_diff(real(vonkarm, r4_kind), ustar_min, neutral,               &
                        & stable_option, new_mo_option, real(rich_crit, r4_kind), &
                        & real(zeta_trans, r4_kind), ni, nj, nk, z_a, u_star_a,   & !miz
                        & b_star_a, k_m_a, k_h_a, ier)
k_m = k_m_a(1,1,1)
k_h = k_h_a(1,1,1)

end subroutine mo_diff_0d_1_r4

!=======================================================================

subroutine mo_diff_0d_n_r4(z, u_star, b_star, k_m, k_h)

real(kind=r4_kind), intent(in),  dimension(:) :: z
real(kind=r4_kind), intent(in)                :: u_star, b_star
real(kind=r4_kind), intent(out), dimension(:) :: k_m, k_h

integer                                            :: ni, nj, nk, ier
integer, parameter                                 :: lkind = r4_kind
real(kind=r4_kind), parameter                 :: ustar_min = 1.0E-10_lkind
real(kind=r4_kind), dimension(1,1,size(z))    :: z_a, k_m_a, k_h_a
real(kind=r4_kind), dimension(1,1)            :: u_star_a, b_star_a

if(.not.module_is_initialized) call error_mesg('mo_diff_0d_n in monin_obukhov_mod', &
     'monin_obukhov_init has not been called', FATAL)

ni = 1; nj = 1; nk = size(z(:))
z_a(1,1,:)    = z(:)
u_star_a(1,1) = u_star
b_star_a(1,1) = b_star
call monin_obukhov_diff(real(vonkarm, r4_kind), ustar_min, neutral,              &
                       & stable_option, new_mo_option, real(rich_crit, r4_kind), &
                       & real(zeta_trans, r4_kind), ni, nj, nk, z_a, u_star_a,   &
                       & b_star_a, k_m_a, k_h_a, ier)
k_m(:) = k_m_a(1,1,:)
k_h(:) = k_h_a(1,1,:)
end subroutine mo_diff_0d_n_r4

!=======================================================================

subroutine stable_mix_2d_r4(rich, mix)

real(kind=r4_kind), intent(in) , dimension(:,:)  :: rich
real(kind=r4_kind), intent(out), dimension(:,:)  :: mix
integer :: n2 !< Size of dimension 2 of mix and rich
integer :: i !< Loop index

n2 = size(mix, 2)

do i=1, n2
  call stable_mix(rich(:, i), mix(:, i))
enddo

end subroutine stable_mix_2d_r4


!=======================================================================

subroutine stable_mix_1d_r4(rich, mix)

real(kind=r4_kind), intent(in) , dimension(:)  :: rich
real(kind=r4_kind), intent(out), dimension(:)  :: mix
integer :: n !< Size of mix and rich
integer :: ierr !< Error code set by monin_obukhov_stable_mix

if (.not.module_is_initialized) call error_mesg('stable_mix in monin_obukhov_mod', &
     'monin_obukhov_init has not been called', FATAL)

n = size(mix)

call monin_obukhov_stable_mix(stable_option, real(rich_crit,r4_kind), &
                              & real(zeta_trans,r4_kind), n, rich, mix, ierr)

end subroutine stable_mix_1d_r4

!=======================================================================

subroutine stable_mix_0d_r4(rich, mix)

real(kind=r4_kind), intent(in)       :: rich
real(kind=r4_kind), intent(out)      :: mix
real(kind=r4_kind), dimension(1)     :: mix_1d !< Representation of mix as a dimension(1) array

call stable_mix([rich], mix_1d)

mix = mix_1d(1)

end subroutine stable_mix_0d_r4
!=======================================================================
# 96 "monin_obukhov/include/monin_obukhov_r4.fh" 2
# 210 "monin_obukhov/monin_obukhov.F90" 2

# 1 "monin_obukhov/include/monin_obukhov_r8.fh" 1
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
















































































# 1 "monin_obukhov/include/monin_obukhov.inc" 1
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


subroutine mo_drag_1d_r8 &
         (pt, pt0, z, z0, zt, zq, speed, drag_m, drag_t, drag_q, &
          u_star, b_star, avail)

real(kind=r8_kind), intent(in)   , dimension(:) :: pt, pt0, z, z0, zt, zq, speed
real(kind=r8_kind), intent(inout), dimension(:) :: drag_m, drag_t, drag_q, u_star, b_star
logical, intent(in), optional, dimension(:)          :: avail

logical                                              :: lavail, avail_dummy(1)
integer                                              :: n, ier

integer, parameter                                   :: max_iter = 20
integer, parameter                                   :: lkind    = r8_kind
real(kind=r8_kind), parameter                   :: error    = 1.0E-04_lkind, &
                                                        zeta_min = 1.0E-06_lkind, &
                                                        small    = 1.0E-04_lkind

! #include "monin_obukhov_interfaces.h"

if(.not.module_is_initialized) call error_mesg('mo_drag_1d in monin_obukhov_mod', &
     'monin_obukhov_init has not been called', FATAL)

n      = size(pt)
lavail = .false.
if(present(avail)) lavail = .true.


if(lavail) then
   if (count(avail) .eq. 0) return
   call monin_obukhov_drag_1d(real(grav, r8_kind), real(vonkarm, r8_kind),    &
        & error, zeta_min, max_iter, real(small, r8_kind), neutral, stable_option, &
        & new_mo_option, real(rich_crit, r8_kind), real(zeta_trans, r8_kind), &!miz
        & real(drag_min_heat, r8_kind), real(drag_min_moist, r8_kind),        &
        & real(drag_min_mom, r8_kind), n, pt, pt0, z, z0, zt,                      &
        & zq, speed, drag_m, drag_t, drag_q, u_star, b_star, lavail, avail, ier)
else
call monin_obukhov_drag_1d(real(grav, r8_kind), real(vonkarm, r8_kind),             &
        & error, zeta_min, max_iter, real(small, r8_kind), neutral, stable_option,       &
        & new_mo_option, real(rich_crit, r8_kind), real(zeta_trans, r8_kind),       &!miz
        & real(drag_min_heat, r8_kind), real(drag_min_moist, r8_kind),              &
        & real(drag_min_mom, r8_kind), n, pt, pt0, z, z0, zt, zq, speed, drag_m, drag_t, &
        & drag_q, u_star, b_star, lavail, avail_dummy, ier)
endif

end subroutine mo_drag_1d_r8


!=======================================================================

subroutine mo_profile_1d_r8(zref, zref_t, z, z0, zt, zq, u_star, b_star, q_star, &
                         del_m, del_t, del_q, avail)

real(kind=r8_kind),    intent(in)                :: zref, zref_t
real(kind=r8_kind),    intent(in) , dimension(:) :: z, z0, zt, zq, u_star, b_star, q_star
real(kind=r8_kind),    intent(out), dimension(:) :: del_m, del_t, del_q
logical, intent(in) , optional, dimension(:)          :: avail

logical                                               :: dummy_avail(1)
integer                                               :: n, ier

! #include "monin_obukhov_interfaces.h"

if(.not.module_is_initialized) call error_mesg('mo_profile_1d in monin_obukhov_mod', &
     'monin_obukhov_init has not been called', FATAL)

n = size(z)
if(present(avail)) then

   if (count(avail) .eq. 0) return

   call monin_obukhov_profile_1d(real(vonkarm, r8_kind), &
        & neutral, stable_option, new_mo_option, real(rich_crit, r8_kind),   &
        & real(zeta_trans, r8_kind), n, zref, zref_t, z, z0, zt, zq, u_star, &
        & b_star, q_star, del_m, del_t, del_q, .true., avail, ier)

else

   call monin_obukhov_profile_1d(real(vonkarm, r8_kind), &
        & neutral, stable_option, new_mo_option, real(rich_crit, r8_kind),   &
        & real(zeta_trans, r8_kind), n, zref, zref_t, z, z0, zt, zq, u_star, &
        & b_star, q_star, del_m, del_t, del_q, .false., dummy_avail, ier)

endif

end subroutine mo_profile_1d_r8

!=======================================================================

subroutine stable_mix_3d_r8(rich, mix)

real(kind=r8_kind), intent(in) , dimension(:,:,:)  :: rich
real(kind=r8_kind), intent(out), dimension(:,:,:)  :: mix
integer :: n2 !< Size of dimension 2 of mix and rich
integer :: n3 !< Size of dimension 3 of mix and rich
integer :: i, j !< Loop indices

n2 = size(mix, 2)
n3 = size(mix, 3)

do j=1, n3
  do i=1, n2
    call stable_mix(rich(:, i, j), mix(:, i, j))
  enddo
enddo

end subroutine stable_mix_3d_r8

!=======================================================================

subroutine mo_diff_2d_n_r8(z, u_star, b_star, k_m, k_h)

real(kind=r8_kind), intent(in),  dimension(:,:,:) :: z
real(kind=r8_kind), intent(in),  dimension(:,:)   :: u_star, b_star
real(kind=r8_kind), intent(out), dimension(:,:,:) :: k_m, k_h

integer                                                :: ni, nj, nk, ier
integer, parameter                                     :: lkind = r8_kind
real(kind=r8_kind), parameter                     :: ustar_min = 1.0E-10_lkind

if(.not.module_is_initialized) call error_mesg('mo_diff_2d_n in monin_obukhov_mod', &
     'monin_obukhov_init has not been called', FATAL)

ni = size(z, 1); nj = size(z, 2); nk = size(z, 3)
call monin_obukhov_diff(real(vonkarm, r8_kind), ustar_min, neutral,                      &
                             & stable_option, new_mo_option, real(rich_crit, r8_kind),   &
                             & real(zeta_trans, r8_kind), ni, nj, nk, z, u_star, b_star, &
                             & k_m, k_h, ier)

end subroutine mo_diff_2d_n_r8

!=======================================================================
! The following routines are used by the public interfaces above
!=======================================================================

subroutine solve_zeta_r8(rich, z, z0, zt, zq, f_m, f_t, f_q, mask)

real(kind=r8_kind), intent(in) , dimension(:) :: rich, z, z0, zt, zq
logical, intent(in) , dimension(:)                 :: mask
real(kind=r8_kind), intent(out), dimension(:) :: f_m, f_t, f_q

integer, parameter                                 :: lkind    = r8_kind
real(kind=r8_kind), parameter                 :: error    = 1.0E-04_lkind
real(kind=r8_kind), parameter                 :: zeta_min = 1.0E-06_lkind
integer, parameter                                 :: max_iter = 20

real(kind=r8_kind)                            :: max_cor
integer                                            :: iter

real(kind=r8_kind), dimension(size(rich(:)))  ::   &
          d_rich, rich_1, correction, corr, z_z0, z_zt, z_zq, &
          ln_z_z0, ln_z_zt, ln_z_zq, zeta,                    &
          phi_m, phi_m_0, phi_t, phi_t_0, rzeta,              &
          zeta_0, zeta_t, zeta_q, df_m, df_t

logical, dimension(size(rich(:)))                  :: mask_1

z_z0    = z/z0
z_zt    = z/zt
z_zq    = z/zq
ln_z_z0 = log(z_z0)
ln_z_zt = log(z_zt)
ln_z_zq = log(z_zq)

corr = 0.0_lkind
mask_1 = mask

! initial guess

where(mask_1)
  zeta = rich*ln_z_z0*ln_z_z0/ln_z_zt
elsewhere
  zeta = 0.0_lkind
end where

where (mask_1 .and. rich >= 0.0_lkind)
  zeta = zeta/(1.0_lkind - rich/real(rich_crit, r8_kind))
end where

iter_loop: do iter = 1, max_iter

  where (mask_1 .and. abs(zeta).lt.zeta_min)
    zeta   = 0.0_lkind
    f_m    = ln_z_z0
    f_t    = ln_z_zt
    f_q    = ln_z_zq
    mask_1 = .false.  ! don't do any more calculations at these pts
  end where

  where (mask_1)
    rzeta  = 1.0_lkind/zeta
    zeta_0 = zeta/z_z0
    zeta_t = zeta/z_zt
    zeta_q = zeta/z_zq
  elsewhere
    zeta_0 = 0.0_lkind
    zeta_t = 0.0_lkind
    zeta_q = 0.0_lkind
  end where

  call mo_derivative_m(phi_m  , zeta  , mask_1)
  call mo_derivative_m(phi_m_0, zeta_0, mask_1)
  call mo_derivative_t(phi_t  , zeta  , mask_1)
  call mo_derivative_t(phi_t_0, zeta_t, mask_1)

  call mo_integral_m(f_m, zeta, zeta_0, ln_z_z0, mask_1)
  call mo_integral_tq(f_t, f_q, zeta, zeta_t, zeta_q, ln_z_zt, ln_z_zq, mask_1)

  where (mask_1)
    df_m       = (phi_m - phi_m_0)*rzeta
    df_t       = (phi_t - phi_t_0)*rzeta
    rich_1     = zeta*f_t/(f_m*f_m)
    d_rich     = rich_1*( rzeta +  df_t/f_t - 2.0_lkind *df_m/f_m)
    correction = (rich - rich_1)/d_rich
    corr       = min(abs(correction),abs(correction/zeta))
      ! the criterion corr < error seems to work ok, but is a bit arbitrary
      !  when zeta is small the tolerance is reduced
  end where

  max_cor = maxval(corr)

  if(max_cor > error) then
    mask_1 = mask_1 .and. (corr > error)
       ! change the mask so computation proceeds only on non-converged points
    where(mask_1)
      zeta = zeta + correction
    end where
    cycle iter_loop
  else
    return
  end if

end do iter_loop

call error_mesg ('solve_zeta in monin_obukhov_mod',  &
                 'surface drag iteration did not converge', FATAL)

end subroutine solve_zeta_r8

!=======================================================================

subroutine mo_derivative_m_r8(phi_m, zeta, mask)

! the differential similarity function for momentum

real(kind=r8_kind), intent(out),  dimension(:) :: phi_m
real(kind=r8_kind), intent(in),   dimension(:) :: zeta
logical,                 intent(in),   dimension(:) :: mask

logical, dimension(size(zeta(:)))                   :: stable, unstable
real(kind=r8_kind), dimension(size(zeta(:)))   :: x
integer, parameter                                  :: lkind = r8_kind

stable   = mask .and. zeta >= 0.0_lkind
unstable = mask .and. zeta <  0.0_lkind

where (unstable)
  x     = (1.0_lkind - 16.0_lkind*zeta  )**(-0.5_lkind)
  phi_m = sqrt(x)  ! phi_m = (1 - 16.0*zeta)**(-0.25)
end where

if(stable_option == 1) then

  where (stable)
    phi_m = 1.0_lkind + zeta*(5.0_lkind + real(b_stab, r8_kind) &
            *zeta)/(1.0_lkind + zeta)
  end where

else if(stable_option == 2) then

  where (stable .and. zeta < real(zeta_trans,r8_kind))
    phi_m = 1.0_lkind + 5.0_lkind*zeta
  end where
  where (stable .and. zeta >= real(zeta_trans,r8_kind))
    phi_m = real(lambda, r8_kind) + real(b_stab, r8_kind)*zeta
  end where

endif

return
end subroutine mo_derivative_m_r8

!=======================================================================

subroutine mo_derivative_t_r8(phi_t, zeta, mask)

! the differential similarity function for buoyancy and tracers

real(kind=r8_kind), intent(out),  dimension(:) :: phi_t
real(kind=r8_kind), intent(in),   dimension(:) :: zeta
logical ,                intent(in),   dimension(:) :: mask

logical, dimension(size(zeta(:)))                   :: stable, unstable
integer, parameter                                  :: lkind = r8_kind

stable   = mask .and. zeta >= 0.0_lkind
unstable = mask .and. zeta <  0.0_lkind

where (unstable)
  phi_t = (1.0_lkind - 16.0_lkind*zeta)**(-0.5_lkind)
end where

if(stable_option == 1) then

  where (stable)
    phi_t = 1.0_lkind + zeta * (5.0_lkind + real(b_stab, r8_kind)&
            * zeta)/(1.0_lkind + zeta)
  end where

else if(stable_option == 2) then

  where (stable .and. zeta < real(zeta_trans,r8_kind))
    phi_t = 1.0_lkind + 5.0_lkind*zeta
  end where
  where (stable .and. zeta >= real(zeta_trans,r8_kind))
    phi_t = real(lambda, r8_kind) + real(b_stab, r8_kind)*zeta
  end where

endif

return
end subroutine mo_derivative_t_r8

!=======================================================================

subroutine mo_integral_tq_r8 (psi_t, psi_q, zeta, zeta_t, zeta_q, &
                           ln_z_zt, ln_z_zq, mask)

! the integral similarity function for moisture and tracers

real(kind=r8_kind), intent(out), dimension(:) :: psi_t, psi_q
real(kind=r8_kind), intent(in),  dimension(:) :: zeta, zeta_t, zeta_q, ln_z_zt, ln_z_zq
logical , intent(in),  dimension(:)                :: mask

real(kind=r8_kind), dimension(size(zeta(:)))  :: x, x_t, x_q

logical, dimension(size(zeta(:)))                  :: stable, unstable, &
                                                      weakly_stable, strongly_stable
integer, parameter                                 :: lkind = r8_kind

stable   = mask .and. zeta >= 0.0_lkind
unstable = mask .and. zeta <  0.0_lkind

where(unstable)

  x     = sqrt(1.0_lkind - 16.0_lkind*zeta)
  x_t   = sqrt(1.0_lkind - 16.0_lkind*zeta_t)
  x_q   = sqrt(1.0_lkind - 16.0_lkind*zeta_q)

  psi_t = ln_z_zt - 2.0_lkind*log( (1.0_lkind + x)/(1.0_lkind + x_t) )
  psi_q = ln_z_zq - 2.0_lkind*log( (1.0_lkind + x)/(1.0_lkind + x_q) )

end where

if( stable_option == 1) then

  where (stable)

    psi_t = ln_z_zt + (5.0_lkind - real(b_stab, r8_kind)) &
            *log((1.0_lkind + zeta)/(1.0_lkind + zeta_t)) &
            + real(b_stab, r8_kind)*(zeta - zeta_t)
    psi_q = ln_z_zq + (5.0_lkind - real(b_stab, r8_kind)) &
            *log((1.0_lkind + zeta)/(1.0_lkind + zeta_q)) &
            + real(b_stab, r8_kind)*(zeta - zeta_q)

  end where

else if (stable_option == 2) then

  weakly_stable   = stable .and. zeta <= real(zeta_trans,r8_kind)
  strongly_stable = stable .and. zeta >  real(zeta_trans,r8_kind)

  where (weakly_stable)
    psi_t = ln_z_zt + 5.0_lkind*(zeta - zeta_t)
    psi_q = ln_z_zq + 5.0_lkind*(zeta - zeta_q)
  end where

  where(strongly_stable)
    x = (real(lambda, r8_kind) - 1.0_lkind)*log(zeta/real(zeta_trans, r8_kind)) + &
         real(b_stab, r8_kind)*(zeta - real(zeta_trans, r8_kind))
  endwhere

  where (strongly_stable .and. zeta_t <= real(zeta_trans,r8_kind))
    psi_t = ln_z_zt + x + 5.0_lkind * (real(zeta_trans, r8_kind) - zeta_t)
  end where

  where (strongly_stable .and. zeta_t > real(zeta_trans,r8_kind))
    psi_t = real(lambda, r8_kind)* ln_z_zt &
            + real(b_stab, r8_kind)*(zeta - zeta_t)
  endwhere

  where (strongly_stable .and. zeta_q <= real(zeta_trans,r8_kind))
    psi_q = ln_z_zq + x + 5.0_lkind &
            *(real(zeta_trans, r8_kind) - zeta_q)
  end where

  where (strongly_stable .and. zeta_q > real(zeta_trans,r8_kind))
    psi_q = real(lambda, r8_kind)*ln_z_zq + real(b_stab, r8_kind) &
            * (zeta - zeta_q)
  endwhere

end if

return
end subroutine mo_integral_tq_r8

!=======================================================================

subroutine mo_integral_m_r8 (psi_m, zeta, zeta_0, ln_z_z0, mask)

!  the integral similarity function for momentum

real(kind=r8_kind), intent(out), dimension(:) :: psi_m
real(kind=r8_kind), intent(in),  dimension(:) :: zeta, zeta_0, ln_z_z0
logical,                 intent(in),  dimension(:) :: mask

real(kind=r8_kind), dimension(size(zeta(:)))  :: x, x_0, x1, x1_0, num, denom, y

logical, dimension(size(zeta(:)))                  :: stable, unstable, &
                                                      weakly_stable, strongly_stable
integer, parameter                                 :: lkind = r8_kind

stable   = mask .and. zeta >= 0.0_lkind
unstable = mask .and. zeta <  0.0_lkind

where(unstable)

  x      = sqrt(1.0_lkind - 16.0_lkind*zeta)
  x_0    = sqrt(1.0_lkind - 16.0_lkind*zeta_0)

  x      = sqrt(x)
  x_0    = sqrt(x_0)

  x1     = 1.0_lkind + x
  x1_0   = 1.0_lkind + x_0

  num    = x1*x1*(1.0_lkind + x*x)
  denom  = x1_0*x1_0*(1.0_lkind + x_0*x_0)
  y      = atan(x) - atan(x_0)
  psi_m  = ln_z_z0 - log(num/denom) + 2.0_lkind*y

end where

if( stable_option == 1) then

  where (stable)
    psi_m = ln_z_z0 + (5.0_lkind - real(b_stab, r8_kind)) &
            *log((1.0_lkind + zeta)/(1.0_lkind + zeta_0)) &
            + real(b_stab, r8_kind)*(zeta - zeta_0)
  end where

else if (stable_option == 2) then

  weakly_stable   = stable .and. zeta <= real(zeta_trans,r8_kind)
  strongly_stable = stable .and. zeta >  real(zeta_trans,r8_kind)

  where (weakly_stable)
    psi_m = ln_z_z0 + 5.0_lkind*(zeta - zeta_0)
  end where

  where(strongly_stable)
    x = (real(lambda, r8_kind) - 1.0_lkind)*log(zeta/real(zeta_trans, r8_kind)) + &
        real(b_stab, r8_kind)*(zeta - real(zeta_trans, r8_kind))
  endwhere

  where (strongly_stable .and. zeta_0 <= real(zeta_trans,r8_kind))
    psi_m = ln_z_z0 + x + 5.0_lkind &
            *(real(zeta_trans, r8_kind) - zeta_0)
  end where
  where (strongly_stable .and. zeta_0 > real(zeta_trans,r8_kind))
    psi_m = real(lambda, r8_kind)*ln_z_z0 + real(b_stab, r8_kind) &
            *(zeta - zeta_0)
  endwhere

end if

return
end subroutine mo_integral_m_r8


!=======================================================================
! The following routines allow the public interfaces to be used
! with different dimensions of the input and output
!
!=======================================================================


subroutine mo_drag_2d_r8 &
    (pt, pt0, z, z0, zt, zq, speed, drag_m, drag_t, drag_q, u_star, b_star)

real(kind=r8_kind), intent(in)   , dimension(:,:) :: z, speed, pt, pt0, z0, zt, zq
real(kind=r8_kind), intent(out)  , dimension(:,:) :: drag_m, drag_t, drag_q
real(kind=r8_kind), intent(inout), dimension(:,:) :: u_star, b_star

integer :: j

do j = 1, size(pt,2)
  call mo_drag (pt(:,j), pt0(:,j), z(:,j), z0(:,j), zt(:,j), zq(:,j), &
                   speed(:,j), drag_m(:,j), drag_t(:,j), drag_q(:,j), &
                   u_star(:,j), b_star(:,j))
end do


return
end subroutine mo_drag_2d_r8

!=======================================================================
subroutine mo_drag_0d_r8 &
    (pt, pt0, z, z0, zt, zq, speed, drag_m, drag_t, drag_q, u_star, b_star)

real(kind=r8_kind), intent(in)    :: z, speed, pt, pt0, z0, zt, zq
real(kind=r8_kind), intent(out)   :: drag_m, drag_t, drag_q, u_star, b_star

real(kind=r8_kind), dimension(1)  :: pt_1, pt0_1, z_1, z0_1, zt_1, zq_1, speed_1, &
                      drag_m_1, drag_t_1, drag_q_1, u_star_1, b_star_1

pt_1   (1) = pt
pt0_1  (1) = pt0
z_1    (1) = z
z0_1   (1) = z0
zt_1   (1) = zt
zq_1   (1) = zq
speed_1(1) = speed

call mo_drag (pt_1, pt0_1, z_1, z0_1, zt_1, zq_1, speed_1, &
                 drag_m_1, drag_t_1, drag_q_1, u_star_1, b_star_1)

drag_m = drag_m_1(1)
drag_t = drag_t_1(1)
drag_q = drag_q_1(1)
u_star = u_star_1(1)
b_star = b_star_1(1)

return
end subroutine mo_drag_0d_r8
!=======================================================================

subroutine mo_profile_2d_r8(zref, zref_t, z, z0, zt, zq, u_star, b_star, q_star, &
                         del_m, del_h, del_q)

real(kind=r8_kind), intent(in)                  :: zref, zref_t
real(kind=r8_kind), intent(in) , dimension(:,:) :: z, z0, zt, zq, u_star, b_star, q_star
real(kind=r8_kind), intent(out), dimension(:,:) :: del_m, del_h, del_q

integer :: j

do j = 1, size(z,2)
  call mo_profile (zref, zref_t, z(:,j), z0(:,j), zt(:,j),         &
                      zq(:,j), u_star(:,j), b_star(:,j), q_star(:,j), &
                      del_m(:,j), del_h (:,j), del_q (:,j))
enddo

return
end subroutine mo_profile_2d_r8

!=======================================================================

subroutine mo_profile_0d_r8(zref, zref_t, z, z0, zt, zq, u_star, b_star, q_star, &
                         del_m, del_h, del_q)

real(kind=r8_kind), intent(in)  :: zref, zref_t
real(kind=r8_kind), intent(in)  :: z, z0, zt, zq, u_star, b_star, q_star
real(kind=r8_kind), intent(out) :: del_m, del_h, del_q

real(kind=r8_kind), dimension(1) :: z_1, z0_1, zt_1, zq_1, u_star_1, b_star_1, q_star_1, &
                      del_m_1, del_h_1, del_q_1

z_1     (1) = z
z0_1    (1) = z0
zt_1    (1) = zt
zq_1    (1) = zq
u_star_1(1) = u_star
b_star_1(1) = b_star
q_star_1(1) = q_star

call mo_profile (zref, zref_t, z_1, z0_1, zt_1, zq_1, &
                    u_star_1, b_star_1, q_star_1,        &
                    del_m_1, del_h_1, del_q_1)

del_m = del_m_1(1)
del_h = del_h_1(1)
del_q = del_q_1(1)


return
end subroutine mo_profile_0d_r8

!=======================================================================

subroutine mo_profile_1d_n_r8(zref, z, z0, zt, zq, u_star, b_star, q_star, &
                         del_m, del_t, del_q, avail)

real(kind=r8_kind),    intent(in),  dimension(:)   :: zref
real(kind=r8_kind),    intent(in) , dimension(:)   :: z, z0, zt, zq, u_star, b_star, q_star
real(kind=r8_kind),    intent(out), dimension(:,:) :: del_m, del_t, del_q
logical, intent(in) , optional,          dimension(:)   :: avail

integer :: k

do k = 1, size(zref(:))
  if(present(avail)) then
    call mo_profile (zref(k), zref(k), z, z0, zt, zq, &
       u_star, b_star, q_star, del_m(:,k), del_t(:,k), del_q(:,k), avail)
  else
      call mo_profile (zref(k), zref(k), z, z0, zt, zq, &
       u_star, b_star, q_star, del_m(:,k), del_t(:,k), del_q(:,k))
  endif
enddo

return
end subroutine mo_profile_1d_n_r8

!=======================================================================

subroutine mo_profile_0d_n_r8(zref, z, z0, zt, zq, u_star, b_star, q_star, &
                         del_m, del_t, del_q)

real(kind=r8_kind),    intent(in),  dimension(:) :: zref
real(kind=r8_kind),    intent(in)                :: z, z0, zt, zq, u_star, b_star, q_star
real(kind=r8_kind),    intent(out), dimension(:) :: del_m, del_t, del_q

integer :: k

do k = 1, size(zref(:))
  call mo_profile (zref(k), zref(k), z, z0, zt, zq, &
       u_star, b_star, q_star, del_m(k), del_t(k), del_q(k))
enddo

return
end subroutine mo_profile_0d_n_r8

!=======================================================================

subroutine mo_profile_2d_n_r8(zref, z, z0, zt, zq, u_star, b_star, q_star, &
                         del_m, del_t, del_q)

real(kind=r8_kind),    intent(in),  dimension(:)     :: zref
real(kind=r8_kind),    intent(in),  dimension(:,:)   :: z, z0, zt, zq, u_star, b_star, q_star
real(kind=r8_kind),    intent(out), dimension(:,:,:) :: del_m, del_t, del_q

integer :: k

do k = 1, size(zref(:))
  call mo_profile (zref(k), zref(k), z, z0, zt, zq, &
       u_star, b_star, q_star, del_m(:,:,k), del_t(:,:,k), del_q(:,:,k))
enddo

return
end subroutine mo_profile_2d_n_r8

!=======================================================================

subroutine mo_diff_2d_1_r8(z, u_star, b_star, k_m, k_h)

real(kind=r8_kind), intent(in),  dimension(:,:)      :: z, u_star, b_star
real(kind=r8_kind), intent(out), dimension(:,:)      :: k_m, k_h

real(kind=r8_kind), dimension(size(z,1),size(z,2),1) :: z_n, k_m_n, k_h_n

z_n(:,:,1) = z

call mo_diff(z_n, u_star, b_star, k_m_n, k_h_n)

k_m = k_m_n(:,:,1)
k_h = k_h_n(:,:,1)

return
end subroutine mo_diff_2d_1_r8


!=======================================================================

subroutine mo_diff_1d_1_r8(z, u_star, b_star, k_m, k_h)

real(kind=r8_kind), intent(in),  dimension(:) :: z, u_star, b_star
real(kind=r8_kind), intent(out), dimension(:) :: k_m, k_h

real(kind=r8_kind), dimension(size(z),1,1)    :: z_n, k_m_n, k_h_n
real(kind=r8_kind), dimension(size(z),1)      :: u_star_n, b_star_n

z_n   (:,1,1) = z
u_star_n(:,1) = u_star
b_star_n(:,1) = b_star

call mo_diff(z_n, u_star_n, b_star_n, k_m_n, k_h_n)

k_m = k_m_n(:,1,1)
k_h = k_h_n(:,1,1)

return
end subroutine mo_diff_1d_1_r8

!=======================================================================

subroutine mo_diff_1d_n_r8(z, u_star, b_star, k_m, k_h)

real(kind=r8_kind), intent(in),  dimension(:,:) :: z
real(kind=r8_kind), intent(in),  dimension(:)   :: u_star, b_star
real(kind=r8_kind), intent(out), dimension(:,:) :: k_m, k_h

real(kind=r8_kind), dimension(size(z,1),1)            :: u_star2, b_star2
real(kind=r8_kind), dimension(size(z,1),1, size(z,2)) :: z2, k_m2, k_h2

integer :: n

do n = 1, size(z,2)
  z2   (:,1,n) = z(:,n)
enddo
u_star2(:,1) = u_star
b_star2(:,1) = b_star

call mo_diff(z2, u_star2, b_star2, k_m2, k_h2)

do n = 1, size(z,2)
  k_m(:,n) = k_m2(:,1,n)
  k_h(:,n) = k_h2(:,1,n)
enddo

return
end subroutine mo_diff_1d_n_r8

!=======================================================================

subroutine mo_diff_0d_1_r8(z, u_star, b_star, k_m, k_h)

real(kind=r8_kind), intent(in)       :: z, u_star, b_star
real(kind=r8_kind), intent(out)      :: k_m, k_h

integer                                   :: ni, nj, nk, ier
integer, parameter                        :: lkind = r8_kind
real(kind=r8_kind), parameter        :: ustar_min = 1.0E-10_lkind
real(kind=r8_kind), dimension(1,1,1) :: z_a, k_m_a, k_h_a
real(kind=r8_kind), dimension(1,1)   :: u_star_a, b_star_a

if(.not.module_is_initialized) call error_mesg('mo_diff_0d_1 in monin_obukhov_mod', &
     'monin_obukhov_init has not been called', FATAL)

ni = 1; nj = 1; nk = 1
z_a(1,1,1)    = z
u_star_a(1,1) = u_star
b_star_a(1,1) = b_star
call monin_obukhov_diff(real(vonkarm, r8_kind), ustar_min, neutral,               &
                        & stable_option, new_mo_option, real(rich_crit, r8_kind), &
                        & real(zeta_trans, r8_kind), ni, nj, nk, z_a, u_star_a,   & !miz
                        & b_star_a, k_m_a, k_h_a, ier)
k_m = k_m_a(1,1,1)
k_h = k_h_a(1,1,1)

end subroutine mo_diff_0d_1_r8

!=======================================================================

subroutine mo_diff_0d_n_r8(z, u_star, b_star, k_m, k_h)

real(kind=r8_kind), intent(in),  dimension(:) :: z
real(kind=r8_kind), intent(in)                :: u_star, b_star
real(kind=r8_kind), intent(out), dimension(:) :: k_m, k_h

integer                                            :: ni, nj, nk, ier
integer, parameter                                 :: lkind = r8_kind
real(kind=r8_kind), parameter                 :: ustar_min = 1.0E-10_lkind
real(kind=r8_kind), dimension(1,1,size(z))    :: z_a, k_m_a, k_h_a
real(kind=r8_kind), dimension(1,1)            :: u_star_a, b_star_a

if(.not.module_is_initialized) call error_mesg('mo_diff_0d_n in monin_obukhov_mod', &
     'monin_obukhov_init has not been called', FATAL)

ni = 1; nj = 1; nk = size(z(:))
z_a(1,1,:)    = z(:)
u_star_a(1,1) = u_star
b_star_a(1,1) = b_star
call monin_obukhov_diff(real(vonkarm, r8_kind), ustar_min, neutral,              &
                       & stable_option, new_mo_option, real(rich_crit, r8_kind), &
                       & real(zeta_trans, r8_kind), ni, nj, nk, z_a, u_star_a,   &
                       & b_star_a, k_m_a, k_h_a, ier)
k_m(:) = k_m_a(1,1,:)
k_h(:) = k_h_a(1,1,:)
end subroutine mo_diff_0d_n_r8

!=======================================================================

subroutine stable_mix_2d_r8(rich, mix)

real(kind=r8_kind), intent(in) , dimension(:,:)  :: rich
real(kind=r8_kind), intent(out), dimension(:,:)  :: mix
integer :: n2 !< Size of dimension 2 of mix and rich
integer :: i !< Loop index

n2 = size(mix, 2)

do i=1, n2
  call stable_mix(rich(:, i), mix(:, i))
enddo

end subroutine stable_mix_2d_r8


!=======================================================================

subroutine stable_mix_1d_r8(rich, mix)

real(kind=r8_kind), intent(in) , dimension(:)  :: rich
real(kind=r8_kind), intent(out), dimension(:)  :: mix
integer :: n !< Size of mix and rich
integer :: ierr !< Error code set by monin_obukhov_stable_mix

if (.not.module_is_initialized) call error_mesg('stable_mix in monin_obukhov_mod', &
     'monin_obukhov_init has not been called', FATAL)

n = size(mix)

call monin_obukhov_stable_mix(stable_option, real(rich_crit,r8_kind), &
                              & real(zeta_trans,r8_kind), n, rich, mix, ierr)

end subroutine stable_mix_1d_r8

!=======================================================================

subroutine stable_mix_0d_r8(rich, mix)

real(kind=r8_kind), intent(in)       :: rich
real(kind=r8_kind), intent(out)      :: mix
real(kind=r8_kind), dimension(1)     :: mix_1d !< Representation of mix as a dimension(1) array

call stable_mix([rich], mix_1d)

mix = mix_1d(1)

end subroutine stable_mix_0d_r8
!=======================================================================
# 100 "monin_obukhov/include/monin_obukhov_r8.fh" 2
# 211 "monin_obukhov/monin_obukhov.F90" 2

end module monin_obukhov_mod
!> @}
! close documentation grouping
