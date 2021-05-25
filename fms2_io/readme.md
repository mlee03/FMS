### Fms_io/mpp_io to fms2_io conversion guide

This guide helps covert fms_io/mpp_io code to fms2_io

### A. FMS2_io Fileobjs

#### 1. FmsNetcdfFile_t
This type provides a thin wrapper over the netCDF4 library, but allows the user to assign a “pelist” to the file. If a pelist is assigned, only the first rank on the list directly interacts with the netCDF library, and performs broadcasts to relay the information to the rest of the ranks on the list.

#### 2. FmsNetcdfDomainFile_t
This type extends the FmsNetcdfFile_t type to support “domain-decomposed” reads/writes on a user-defined mpp_domains two-dimensional lon-lat or cubed-sphere grid.

#### 3. FmsNetcdfUnstructuredDomainFile_t
This type extends the FmsNetcdfFile_t type to support “domain-decomposed” reads/writes on a user defined mpp_domains unstructured grid. 

### B. Writing Restarts

#### 1. Domain Decomposed Restarts
**FMS_io**
```F90

use fms_io_mod,         only: restart_file_type, register_restart_field, save_restart
use mpp_domains_mod,    only: domain

integer                 :: id_restart     !< Id for restart variable
type(restart_file_type) :: fileobj        !< Fms_io fileobj
real, dimension(:,:)    :: variable_data  !< Variable data in the compute or data domain
type (domain2d)         :: domain         !< 2d mpp domain

id_restart = register_restart_field(fileobj, "filename", "variable_name", variable_data, domain=domain)
call save_restart(Atm_restart)
```

Filename mangling:
- If the domain is a cubesphere and your io_layout is (/1,1/) or if your io_domain is not defined, it will create restart files: "RESTART/filename.tileX.res.nc"
- If the domain is a cubesphere and your io_layout is (/X,Y/), it write `X*Y` restart files: "RESTART/filename.tileX.res.nc.XXXX"
- If the domain is not a cubephere or it is only 1 tile, "tileX" will not be added

Metadata:
- FMS_io wrote dimensions:
- FMS_io wrote the dimensions as variables as well 
- FMS_io wrote variable attribute: "longname = {same as variable" and "units = "none"} by default

**FMS2_io**
```F90

use fms2_io_mod,        only: FmsNetcdfDomainFile_t, register_restart_field, register_axis, unlimited
use fms2_io_mod,        only: open_file, close_file, write_restart

type(FmsNetcdfDomainFile_t) ::  fileobj

atm_restart_exist = open_file(fileobj, "filename", "ovewrite", domain, is_restart=.true.)
if (atm_restart_exist) then
  call register_axis(Til_restart, dim_names(1), "x")
  call register_axis(Til_restart, dim_names(2), "y")
  call register_axis(Til_restart, dim_names(3), unlimited)
  
  call register_restart_field(Til_restart, 'lprec', Atmos%lprec, dim_names)
  call write_restart
  call close_file
endif
```
- `open_file`: 
  -  Mangles the filename in the same manner as fms_io
  -  Opens the netcdf file to write (a `nf90_create` call)
  -  Set ups the pelist for each io_domain
  -  `is_restart` indicates that this is a restart file, so it adds ".res" to the filename and it allows user to use the `read_restart` and `write_restart` functionality
- `register_axis` 
  - writes the dimension metadata in the netcdf file (a `nf90_def_dim` call)
  - The "x" and "y" indicate that that dimension is domain decomposed in x/y
  - The "unlimited" indicate that the dimension is unlimited dimension
  - The "" indicate that this is a normal dimension of length "". 
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

Other helpful information:
- By default, fms2_io writes the file as `nf90_64bit_offset`, the user can change the netcdf file type by `nc_format="64bit", "classic", or "netcdf4"` to the `open_file` call
- If the user wishes to not add ".res" to filename, the user can add `dont_add_res_to_filename=.true.` to the `open_file` call
- Variable attributes can be written by calling `register_variable_attribute` scalar and 1d real and integers (r32 and 64 bit) and string values are supported
```F90
call register_variable_attribute(fileobj, "varname", "attribute_name", value)
```
- Global attributes can be written by calling `register_global_attribute` scalar and 1d real and integers (32 and 64 bit) and scalar string values are supported 
```F90
call register_global_attribute(fileobj, "global_attribute_name", value)
```

#### 3. Non-domain Decomposed Restarts
Restart files without domain decomposed variables can be written using the `` fileobj, similar to the fms2_io. 

With the *optional* pelist argument, only the first rank interacts with the file (opens, writes, closes) and broadcasts the information to the rest of the ranks on the list. 

### C. Reading Restarts

### D. Reading non-restarts

### E. Coupler type restarts

### F. Boundary conditions restarts
