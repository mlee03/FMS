# 1 "drifters/drifters.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "drifters/drifters.F90"
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


# 1 "drifters/fms_switches.h" 1
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


!DEC$ MESSAGE:'Compiling in serial mode'
# 20 "drifters/drifters.F90" 2


!> @defgroup drifters_mod drifters_mod
!! @ingroup drifters
!! @{
!! @brief <TT>Drifters_mod</TT>is a module designed to advect a set of particles, in parallel or
!!   sequentially, given an prescribed velocity field.
!! @author Alexander Pletzer
!!
!> Drifters are idealized point particles with positions that evolve in time according
!! to a prescribed velocity field, starting from some initial conditions. Drifters have
!! no mass, no energy, no size, and no friction and therefore have no impact on the
!! dynamics of the underlying system. The only feature that distinguishes a drifter
!! from another is its trajectory. This makes drifters ideal for tracking pollution
!! clouds and probing fields (e.g. temperature, salinity) along ocean currents, to name
!! a few applications.
!! Drifters can mimic real experiments such as the Argo floats
!! http://www.metoffice.com/research/ocean/argo/ukfloats.html.
!!
!! When run in parallel, on a 2d decomposed domain, <TT>drifters_mod</TT> will handle all the
!! bookkeeping and communication transparently for the user. This involves adding/removing
!! drifters as they enter/leave a processor element (PE) domain. Note that the number of drifters
!! can vary greatly both between PE domains and within a PE domain in the course of a simulation; the drifters'
!! module will also manage dynamically the memory for the user.
!!
!! There are a number of basic assumptions which could make the drifters' module
!! ill-suited for some tasks. First and foremost, it is assumed that the motion of
!! drifters is not erratic but follows deterministic trajectories. Furthermore,
!! drifters should not cross both compute and data domain boundaries within less
!! than a time step. This limitation is imposed by the Runge-Kutta integration
!! scheme, which must be able to complete, within a time step, a trajectory
!! calculation that starts inside the compute domain and ends inside the data domain. Therefore, the drifters,
!! as they are presently modelled, are unlikely to work for very fast objects.
!! This constraint also puts a upper limit to the domain decomposition, although
!! it can often be remedied by increasing the number of ghost nodes.
!!
!! Another fundamental assumption is that the (e.g. velocity) fields are structured,
!! on a per PE domain basis. There is no support for locally nested or unstrucured
!! meshes. Meshes need not be smooth and continuous across PE domains, however.

module drifters_mod
# 943 "drifters/drifters.F90"
end module drifters_mod
!> @}
! close documentation grouping
