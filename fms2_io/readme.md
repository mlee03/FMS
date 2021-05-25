### Fms_io/mpp_io to fms2_io conversion guide

This guide helps covert fms_io/mpp_io code to fms2_io

### A. FMS2_io Fileobjs

#### 1. FmsNetcdfFile_t
This type provides a thin wrapper over the netCDF4 library, but allows the user to assign a “pelist” to the file. If a pelist is assigned, only the first rank on the list directly interacts with the netCDF library, and performs broadcasts to relay the information to the rest of the ranks on the list.

#### 2. FmsNetcdfDomainFile_t
This type extends the FmsNetcdfFile_t type to support “domain-decomposed” reads/writes on a user-defined mpp_domains two-dimensional lon-lat or cubed-sphere grid. This fileobj can write/read non-domain-decomposed variables as well. 

#### 3. FmsNetcdfUnstructuredDomainFile_t
This type extends the FmsNetcdfFile_t type to support “domain-decomposed” reads/writes on a user defined mpp_domains unstructured grid. This fileobj can write/read non-domain-decomposed variables as well.

### B. Writing Restarts
 
#### 1. Domain Decomposed Restarts

FMS mangles the filename for domain decomposed restarts in the following way:
- If the domain is a cubesphere and the io_layout is (/1,1/), it will create restart files: "RESTART/filename.tileX.res.nc"
- If the domain is a cubesphere and the io_layout is (/X,Y/), it write `X*Y` restart files: "RESTART/filename.tileX.res.nc.XXXX"
- If the domain is not a cubephere or it is only 1 tile, "tileX" will not be added

**This is how a restart file can be written with FMS_io:**
```F90

use fms_io_mod,         only: restart_file_type, register_restart_field, save_restart
use mpp_domains_mod,    only: domain

integer                 :: id_restart     !< Id for restart variable
type(restart_file_type) :: fileobj        !< Fms_io fileobj
real, dimension(:,:,:)  :: variable_data  !< Variable data in the compute or data domain
type (domain2d)         :: domain         !< 2d mpp domain

id_restart = register_restart_field(fileobj, "filename", "variable_name", variable_data, domain=domain)
call save_restart(Atm_restart)
```
Metadata:
- FMS_io wrote named the dimensions: "xaxis_1", "yaxis_1", zaxis_1" and "Time"
- FMS_io wrote the dimensions as variables as well 
- FMS_io wrote variable attribute: "longname = {same as variable} and "units = {"none"} to all variables by default

**This is how a restart file can be written with FMS2-io:**
```F90

use fms2_io_mod,        only: FmsNetcdfDomainFile_t, register_restart_field, register_axis, unlimited
use fms2_io_mod,        only: open_file, close_file, write_restart
use mpp_domains_mod,    only: domain, center

type(FmsNetcdfDomainFile_t) :: fileobj        !< Fms2_io domain decomposed fileobj
real, dimension(:,:,:,:)    :: variable_data  !< Variable data in the compute or data domain
type (domain2d)             :: domain         !< 2d mpp domain
character(len=8),           :: dim_names(4)   !< Array of dimension names

dim_names(1) = "xaxis_1"
dim_names(2) = "yaxis_1"
dim_names(3) = "zaxis_1"
dim_names(4) = "Time"

if (open_file(fileobj, "filename", "ovewrite", domain, is_restart=.true.)) then
  call register_axis(fileobj, dim_names(1), "x", position=center)
  call register_axis(fileobj, dim_names(2), "y", position=center)
  call register_axis(fileobj, dim_names(3), dimsize)
  call register_axis(fileobj, dim_names(4), unlimited)
  
  call register_restart_field(fileobj, 'variable_name', variable_data, dim_names)
  call write_restart(fileobj)
  call close_file(fileobj)
endif
```
- `open_file`: 
  -  a logical function, outputs .true. if the file was opened successfully and .false. if it failed. 
  -  Mangles the filename in the same manner as fms_io
  -  Opens the netcdf file to write (a `nf90_create` call)
  -  Set ups the pelist for each io_domain
  -  `is_restart` indicates that this is a restart file, so it adds ".res" to the filename and it allows user to use the `read_restart` and `write_restart` functionality
- `register_axis` 
  - writes the dimension metadata in the netcdf file (a `nf90_def_dim` call)
  - The "x" and "y" indicate that that dimension is domain decomposed in x/y
  - The `position=center` indicates position of the axis (this is the default). The other acceptable values are `position=east` for "x" and`position=north` for "y", in this cases the axis will have an extra point. 
  - The "unlimited" indicates that the dimension is unlimited `nf90_unlimited`
  - The integer "dimsize" indicate that this is a normal dimension of length equal to dimsize 
- `register_restart_field` 
  - Writes the variable metadata to the file (a `nf90_def_var` call)
  - Saves the data as pointers, which will be written to the netcdf file later
- `write_restart`
  - Loops to the restart variables registered
  - Calculates and writes a global checksum for each variables
  - Writes the data to the file (a `nf90_put_var` call)
- `close_file` 
  - Cleans up the fileobj
  - Closes the netcdf file (a `nf90_close` call)

#### 3. Non-domain Decomposed Restarts
Restart files without domain decomposed variables can be written using the `FmsNetcdfFile_t` fileobj. 

```F90

use fms2_io_mod,        only: FmsNetcdfFile_t, register_restart_field, register_axis, unlimited
use fms2_io_mod,        only: open_file, close_file, write_restart
use mpp_mod,            only: mpp_npes, mpp_get_current_pelist

type(FmsNetcdfFile_t)       :: fileobj        !< Fms2_io fileobj
real, dimension(:,:,:,:)    :: variable_data  !< Variable data
character(len=8),           :: dim_names(4)   !< Array of dimension names
integer, allocatable,       :: pes(:)         !< Array of the pes in the current pelist

!< Get the current pelist
allocate(pes(mpp_npes()))
call mpp_get_current_pelist(pes)

dim_names(1) = "xaxis_1"
dim_names(2) = "yaxis_1"
dim_names(3) = "zaxis_1"
dim_names(4) = "Time"

if (open_file(fileobj, "filename", "ovewrite", pelist=pes, is_restart=.true.)) then
  call register_axis(fileobj, dim_names(1), 96)
  call register_axis(fileobj, dim_names(2), 96)
  call register_axis(fileobj, dim_names(3), dimsize)
  call register_axis(fileobj, dim_names(4), unlimited)
  
  call register_restart_field(fileobj, 'variable_name', variable_data, dim_names)
  call write_restart(fileobj)
  call close_file(fileobj)
endif
```
Difference from domain decomposed restarts:
- ".res" is the only thing that is added to the filename
- "x" and "y" cannot be added to the `register_axis` calls
- With the *optional* pelist argument, only the first rank interacts with the file (opens, writes, closes) and broadcasts the information to the rest of the ranks on the list. 

#### 4. Other helpful information:
- By default, fms2_io writes the file as `nf90_64bit_offset`, the user can change the netcdf file type by `nc_format="64bit", "classic", or "netcdf4"` to the `open_file` call or by adding `netcdf_default_format="64bit", "classic", or "netcdf4"` to the fms2_io_nml
- If the user wishes to not add ".res" to filename, the user can add `dont_add_res_to_filename=.true.` to the `open_file` call
- **Variable attributes** can be written by calling `register_variable_attribute` scalar and 1d real and integers (32 and 64 bit) and string values are supported
```F90
call register_variable_attribute(fileobj, "varname", "attribute_name", value)
```
- **Global attributes** can be written by calling `register_global_attribute` scalar and 1d real and integers (32 and 64 bit) and scalar string values are supported 
```F90
call register_global_attribute(fileobj, "global_attribute_name", value)
```
### C. Reading Restarts
The restarts can be read the same way as the writes. The only difference is that "read" is used in the `open_file` call and `read_restart` is used instead of `write restart`. 

#### Other helpful information:

### D. Reading non-restarts

### E. Coupler type restarts

### F. Boundary conditions restarts
