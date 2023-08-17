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

!***********************************************************************
!* Requires namelist but it can be empty. Requires a diag_table such as
!* below, requires the input files: (aerosol.climatology.nc and
!* o3.climatology.nc)

!* diag_table
!* test_diag_manager_01
!* 1 3 1 0 0 0

!* #output files
!*  "diag_test_01",  1, "days", 1, "days", "time"

!* #output variables
!* "test_diag_manager_mod", "dat1", "dat1", "diag_test_01",  "all", .false., "none", 2
!***********************************************************************

program test_interpolator

end program test_interpolator
