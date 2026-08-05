# 1 "constants/fmsconstants.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "constants/fmsconstants.F90"
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
!> @defgroup fmsconstants FMSConstants
!! @ingroup constants
!! @{
!! @brief Defines useful constants for Earth. Constants are defined as real
!!
!!    FMSconstants have been declared as REAL(kind=sizeof(rvar)), PARAMETER.
!!
!!    The value of a constant defined and used from here cannot be changed
!!    in a users program. New constants can be defined in terms of values
!!    from the FMSconstants module and their includes using a parameter
!!    statement.<br><br>
!!
!!    The currently support contant systems are:
!!       GFDL constants (gfdl_constants.fh)
!!       GEOS constants (geos_constants.fh)
!!       GFS  constants (gfs_constants.fh)
!!       <br><br>
!!
!!    The name given to a particular constant may be changed.<br><br>
!!
!!    Constants can only be used on the right side on an assignment statement
!!    (their value can not be reassigned).
!!
!!    Example:
!!
!! @verbatim
!!    use FMSConstants, only:  TFREEZE, grav_new => GRAV
!!    real, parameter :: grav_inv = 1.0 / grav_new
!!    tempc(:,:,:) = tempk(:,:,:) - TFREEZE
!!    geopotential(:,:) = height(:,:) * grav_new
!! @endverbatim

module FMSconstants

  use platform_mod, only: r4_kind, r8_kind

  !--- default scoping
  implicit none

  !--- needed with implicit none
  real :: dum  !< dummy real variable



!--- set a default for the FMSConstants




!--- perform error checking and include the correct system of constants



# 1 "constants/gfdl_constants.fh" 1
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

character(len=18), public, parameter :: constants_version = 'FMSConstants: GFDL'

!--- temporary definition for backwards compatibility
real(kind=sizeof(dum)), public, parameter :: small_fac = 1._r8_kind

!--- Spherical coordinate conversion constants
real(kind=r8_kind), public, parameter :: PI_8 = 3.14159265358979323846_r8_kind  !< Ratio of circle circumference to
                                                                                !! diameter [N/A]
real(kind=sizeof(dum)),   public, parameter :: PI   = PI_8                            !< Ratio of circle circumference to
                                                                                !! diameter [N/A]
real(kind=sizeof(dum)),   public, parameter :: RAD_TO_DEG  = 180._r8_kind/PI_8        !< Degrees per radian [deg/rad]
real(kind=sizeof(dum)),   public, parameter :: DEG_TO_RAD  = PI_8/180._r8_kind        !< Radians per degree [rad/deg]
real(kind=sizeof(dum)),   public, parameter :: RADIAN      = RAD_TO_DEG               !< Equal to RAD_TO_DEG for backward
                                                                                !! compatability. [rad/deg]

!--- Earth physical constants
real(kind=sizeof(dum)), public, parameter :: RADIUS             = 6371.0E+3_r8_kind  !< Radius of the Earth [m]
real(kind=sizeof(dum)), public, parameter :: OMEGA              = 7.292E-5_r8_kind   !< Rotation rate of the Earth [1/s]
real(kind=sizeof(dum)), public, parameter :: GRAV               = 9.80_r8_kind       !< Acceleration due to gravity [m/s^2]
real(kind=sizeof(dum)), public, parameter :: SECONDS_PER_DAY    = 86400._r8_kind     !< Seconds in a day [s]
real(kind=sizeof(dum)), public, parameter :: SECONDS_PER_HOUR   =  3600._r8_kind     !< Seconds in an hour [s]
real(kind=sizeof(dum)), public, parameter :: SECONDS_PER_MINUTE =    60._r8_kind     !< Seconds in a minute [s]

!--- Various gas constants
real(kind=sizeof(dum)), public, parameter :: RDGAS    = 287.04_r8_kind            !< Gas constant for dry air [J/kg/deg]
real(kind=sizeof(dum)), public, parameter :: RVGAS    = 461.50_r8_kind            !< Gas constant for water vapor
                                                                            !! [J/kg/deg]
real(kind=sizeof(dum)), public, parameter :: HLV      = 2.500E6_r8_kind           !< Latent heat of evaporation [J/kg]
real(kind=sizeof(dum)), public, parameter :: HLF      = 3.34E5_r8_kind            !< Latent heat of fusion [J/kg]
real(kind=sizeof(dum)), public, parameter :: HLS      = HLV + HLF                 !< Latent heat of sublimation [J/kg]
real(kind=sizeof(dum)), public, parameter :: KAPPA    = 2.0_r8_kind/7.0_r8_kind   !< RDGAS / CP_AIR [dimensionless]
real(kind=sizeof(dum)), public, parameter :: CP_AIR   = RDGAS/KAPPA               !< Specific heat capacity of dry air
                                                                            !! at constant pressure [J/kg/deg]
real(kind=sizeof(dum)), public, parameter :: CP_VAPOR = 4.0_r8_kind*RVGAS         !< Specific heat capacity of water vapor
                                                                            !! at constant pressure [J/kg/deg]
real(kind=sizeof(dum)), public, parameter :: CP_OCEAN = 3989.24495292815_r8_kind  !< Specific heat capacity taken
                                                                            !! from McDougall (2002)
                                                                            !! "Potential Enthalpy ..." [J/kg/deg]
real(kind=sizeof(dum)), public, parameter :: DENS_H2O = 1000._r8_kind             !< Density of liquid water [kg/m^3]
real(kind=sizeof(dum)), public, parameter :: RHOAIR   = 1.292269_r8_kind          !< Reference atmospheric density [kg/m^3]
real(kind=sizeof(dum)), public, parameter :: RHO0     = 1.035E3_r8_kind           !< Average density of sea water [kg/m^3]
real(kind=sizeof(dum)), public, parameter :: RHO0R    = 1.0_r8_kind/RHO0          !< Reciprocal of average density of
                                                                            !! sea water [m^3/kg]
real(kind=sizeof(dum)), public, parameter :: RHO_CP   = RHO0*CP_OCEAN             !< (kg/m^3)*(cal/kg/deg C)(joules/cal) =
                                                                            !! (joules/m^3/deg C) [J/m^3/deg]
real(kind=sizeof(dum)), public, parameter :: O2MIXRAT = 2.0953E-01_r8_kind        !< Mixing ratio of molecular oxygen
                                                                            !! in air [dimensionless]
real(kind=sizeof(dum)), public, parameter :: WTMAIR   = 2.896440E+01_r8_kind      !< Molecular weight of air [AMU]
real(kind=sizeof(dum)), public, parameter :: WTMH2O   = WTMAIR*(RDGAS/RVGAS)      !< Molecular weight of water [AMU]
real(kind=sizeof(dum)), public, parameter :: WTMOZONE =  47.99820_r8_kind         !< Molecular weight of ozone [AMU]
real(kind=sizeof(dum)), public, parameter :: WTMC     =  12.00000_r8_kind         !< Molecular weight of carbon [AMU]
real(kind=sizeof(dum)), public, parameter :: WTMCO2   =  44.00995_r8_kind         !< Molecular weight of carbon dioxide
                                                                            !! [AMU]
real(kind=sizeof(dum)), public, parameter :: WTMCH4   =  16.0425_r8_kind          !< Molecular weight of methane [AMU]
real(kind=sizeof(dum)), public, parameter :: WTMO2    =  31.9988_r8_kind          !< Molecular weight of molecular
                                                                            !! oxygen [AMU]
real(kind=sizeof(dum)), public, parameter :: WTMCFC11 = 137.3681_r8_kind          !< Molecular weight of CFC-11
                                                                            !! (CCl3F) [AMU]
real(kind=sizeof(dum)), public, parameter :: WTMCFC12 = 120.9135_r8_kind          !< Molecular weight of CFC-21
                                                                            !! (CCl2F2) [AMU]
real(kind=sizeof(dum)), public, parameter :: WTMN     =  14.0067_r8_kind          !< Molecular weight of Nitrogen [AMU]
real(kind=sizeof(dum)), public, parameter :: DIFFAC   = 1.660_r8_kind             !< Diffusivity factor [dimensionless]
real(kind=sizeof(dum)), public, parameter :: ES0      = 1.0_r8_kind               !< Humidity factor [dimensionless]
                                                                            !! Controls the humidity content of
                                                                            !! the atmosphere through
                                                                            !! the Saturation Vapour Pressure
                                                                            !! expression when using DO_SIMPLE

!--- Pressure and Temperature constants
real(kind=sizeof(dum)), public, parameter :: PSTD     = 1.013250E+06_r8_kind      !< Mean sea level pressure [dynes/cm^2]
real(kind=sizeof(dum)), public, parameter :: PSTD_MKS = 101325.0_r8_kind          !< Mean sea level pressure [N/m^2]
real(kind=sizeof(dum)), public, parameter :: KELVIN   = 273.15_r8_kind            !< Degrees Kelvin at zero Celsius [K]
real(kind=sizeof(dum)), public, parameter :: TFREEZE  = 273.16_r8_kind            !< Freezing temperature of fresh water [K]
real(kind=sizeof(dum)), public, parameter :: C2DBARS  = 1.E-4_r8_kind             !< Converts rho*g*z (in mks) to dbars:
                                                                            !! 1dbar = 10^4 (kg/m^3)(m/s^2)m [dbars]

!--- Named constants
real(kind=sizeof(dum)), public, parameter :: STEFAN   = 5.6734E-8_r8_kind         !< Stefan-Boltzmann constant [W/m^2/deg^4]
real(kind=sizeof(dum)), public, parameter :: AVOGNO   = 6.023000E+23_r8_kind      !< Avogadro's number [atoms/mole]
real(kind=sizeof(dum)), public, parameter :: VONKARM  = 0.40_r8_kind              !< Von Karman constant [dimensionless]

!--- Miscellaneous constants
real(kind=sizeof(dum)), public, parameter :: ALOGMIN    = -50.0_r8_kind        !< Minimum value allowed as argument
                                                                         !! to log function [N/A]
real(kind=sizeof(dum)), public, parameter :: EPSLN      = 1.0E-40_r8_kind      !< A small number to prevent
                                                                         !! divide by zero exceptions [N/A]
real(kind=sizeof(dum)), public, parameter :: RADCON     = ((1.0E+02*GRAV)/(1.0D+04*CP_AIR))*SECONDS_PER_DAY !< Factor to
                                                                        !! convert flux divergence
                                                                        !! to heating rate in degrees per day
                                                                        !! [deg sec/(cm day)]
real(kind=sizeof(dum)), public, parameter :: RADCON_MKS = (GRAV/CP_AIR)*SECONDS_PER_DAY !< Factor to
                                                                        !! convert flux divergence
                                                                        !! to heating rate in degrees per day
                                                                        !! [deg sec/(m day)]
# 71 "constants/fmsconstants.F90" 2
# 80 "constants/fmsconstants.F90"

  !--- public interfaces
  public :: FMSConstants_init

  contains

    !> @brief FMSconstants init routine
    subroutine FMSconstants_init
      use mpp_mod, only: stdlog
      integer :: logunit
      logunit = stdlog()

      write (logunit,'(/,80("="),/(a))') trim(constants_version)

    end subroutine FMSconstants_init

end module FMSconstants
!> @}
! close documentation grouping
