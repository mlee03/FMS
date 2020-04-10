cmake_minimum_required (VERSION 2.8)

set(FC  $ENV{FC})
set(CC  $ENV{CC})
set(CXX $ENV{CXX})

#: build libFMS with UFS settings
project( libFMS C CXX Fortran )


#: get MPI
find_package( MPI REQUIRED )

#: set NETCDF
set(NETCDF_INC_DIR $ENV{NETCDF}/include)
set(NETCDF_LIBDIR $ENV{NETCDF}/lib)
set(NETCDF_LIBS -L$ENV{NETCDF}/lib -lnetcdff -lnetcdf)

message("****************" )
message( "${NETCDF}" )
message( "${NETCDF_INC_DIR}" )
message("****************" )
