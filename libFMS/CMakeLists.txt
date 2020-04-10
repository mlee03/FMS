cmake_minimum_required (VERSION 2.8)

set(FC  $ENV{FC})
set(CC  $ENV{CC})
set(CXX $ENV{CXX})

#: build libFMS with UFS settings
project( libFMS C CXX Fortran )

#: copied from UFS CMakeLists.txt
#: env variables exported from ./make_libfms.sh
message( "******************")
foreach(env_var IN ITEMS CMAKE_C_COMPILER CMAKE_CXX_COMPILER CMAKE_Fortran_COMPILER CMAKE_Platform)
  message( "$ENV{${env_var}}" )
  if ( NOT DEFINED ENV{${env_var}})
    message(FATAL_ERROR "${env_var} is not defined")
  endif()
endforeach()
message( "******************")

#: copied from UFS CMakeLists.txt
set(CMAKE_C_COMPILER       $ENV{CMAKE_C_COMPILER})
set(CMAKE_CXX_COMPILER     $ENV{CMAKE_CXX_COMPILER})
set(CMAKE_Fortran_COMPILER $ENV{CMAKE_Fortran_COMPILER})
set(CMAKE_Platform         $ENV{CMAKE_Platform})

#: FMS directories
set( UFS_DIR $ENV{UFS_DIR} )
set( CURR_DIR $ENV{CURR_DIR} )
set( FMS_SRC ${CURR_DIR}/src )


#: directory for Fortran modules
message( "******************" )
message( "current directry ${CURR_DIR}" )
message( "UFS directory    ${UFS_DIR}"  )
message( "FMS directory    ${FMS_SRC}"  )
message( "******************" ) 

#: setup the UFS compile environments  
set ( CMAKE_UFS_DIR ${UFS_DIR}/cmake )

#: get MPI
find_package( MPI REQUIRED )

#: set compiler options and netcdf flags.  will complain about other libraries not being found
#: this will also calls Intel.cmake
include( ${CMAKE_UFS_DIR}/configure_${CMAKE_Platform}.cmake )


#: 
add_definitions(-Duse_libMPI)
add_definitions(-Duse_netCDF)
add_definitions(-DSPMD)
add_definitions(-Duse_LARGEFILE)
add_definitions(-DINTERNAL_FILE_NML)

add_library(
    fms

    ${FMS_SRC}/amip_interp/amip_interp.F90
    ${FMS_SRC}/astronomy/astronomy.F90
    ${FMS_SRC}/axis_utils/axis_utils.F90
    ${FMS_SRC}/axis_utils/axis_utils2.F90
    ${FMS_SRC}/block_control/block_control.F90
    ${FMS_SRC}/column_diagnostics/column_diagnostics.F90
    ${FMS_SRC}/constants/constants.F90
    ${FMS_SRC}/coupler/atmos_ocean_fluxes.F90
    ${FMS_SRC}/coupler/coupler_types.F90
    ${FMS_SRC}/coupler/ensemble_manager.F90
    ${FMS_SRC}/data_override/data_override.F90
    ${FMS_SRC}/diag_integral/diag_integral.F90
    ${FMS_SRC}/diag_manager/diag_axis.F90
    ${FMS_SRC}/diag_manager/diag_data.F90
    ${FMS_SRC}/diag_manager/diag_grid.F90
    ${FMS_SRC}/diag_manager/diag_manager.F90
    ${FMS_SRC}/diag_manager/diag_output.F90
    ${FMS_SRC}/diag_manager/diag_table.F90
    ${FMS_SRC}/diag_manager/diag_util.F90
    ${FMS_SRC}/drifters/cloud_interpolator.F90
    ${FMS_SRC}/drifters/drifters_comm.F90
    ${FMS_SRC}/drifters/drifters_core.F90
    ${FMS_SRC}/drifters/drifters.F90
    ${FMS_SRC}/drifters/drifters_input.F90
    ${FMS_SRC}/drifters/drifters_io.F90
    ${FMS_SRC}/drifters/quicksort.F90
    ${FMS_SRC}/exchange/stock_constants.F90
    ${FMS_SRC}/exchange/xgrid.F90
    ${FMS_SRC}/field_manager/field_manager.F90
    ${FMS_SRC}/field_manager/fm_util.F90
    ${FMS_SRC}/fms/fms.F90
    ${FMS_SRC}/fms/fms_io.F90
    ${FMS_SRC}/fms2_io/blackboxio.F90
    ${FMS_SRC}/fms2_io/fms2_io.F90
    ${FMS_SRC}/fms2_io/fms_io_utils.F90
    ${FMS_SRC}/fms2_io/fms_netcdf_domain_io.F90
    ${FMS_SRC}/fms2_io/fms_netcdf_unstructured_domain_io.F90
    ${FMS_SRC}/fms2_io/netcdf_io.F90
    ${FMS_SRC}/horiz_interp/horiz_interp_bicubic.F90
    ${FMS_SRC}/horiz_interp/horiz_interp_bilinear.F90
    ${FMS_SRC}/horiz_interp/horiz_interp_conserve.F90
    ${FMS_SRC}/horiz_interp/horiz_interp.F90
    ${FMS_SRC}/horiz_interp/horiz_interp_spherical.F90
    ${FMS_SRC}/horiz_interp/horiz_interp_type.F90
    ${FMS_SRC}/interpolator/interpolator.F90
    ${FMS_SRC}/memutils/memutils.F90
    ${FMS_SRC}/mosaic/gradient.F90
    ${FMS_SRC}/mosaic/grid.F90
    ${FMS_SRC}/mosaic/mosaic.F90
    ${FMS_SRC}/mosaic/mosaic2.F90
    ${FMS_SRC}/mpp/mpp_data.F90
    ${FMS_SRC}/mpp/mpp_domains.F90
    ${FMS_SRC}/mpp/mpp_efp.F90
    ${FMS_SRC}/mpp/mpp.F90
    ${FMS_SRC}/mpp/mpp_io.F90
    ${FMS_SRC}/mpp/mpp_memutils.F90
    ${FMS_SRC}/mpp/mpp_parameter.F90
    ${FMS_SRC}/mpp/mpp_pset.F90
    ${FMS_SRC}/mpp/mpp_utilities.F90
    ${FMS_SRC}/platform/platform.F90
    ${FMS_SRC}/random_numbers/MersenneTwister.F90
    ${FMS_SRC}/random_numbers/random_numbers.F90
    ${FMS_SRC}/sat_vapor_pres/sat_vapor_pres.F90
    ${FMS_SRC}/sat_vapor_pres/sat_vapor_pres_k.F90
    ${FMS_SRC}/time_interp/time_interp_external.F90
    ${FMS_SRC}/time_interp/time_interp_external2.F90
    ${FMS_SRC}/time_interp/time_interp.F90
    ${FMS_SRC}/time_manager/get_cal_time.F90
    ${FMS_SRC}/time_manager/time_manager.F90
    ${FMS_SRC}/topography/gaussian_topog.F90
    ${FMS_SRC}/topography/topography.F90
    ${FMS_SRC}/tracer_manager/tracer_manager.F90
    ${FMS_SRC}/tridiagonal/tridiagonal.F90
    ${FMS_SRC}/mosaic/create_xgrid.c
    ${FMS_SRC}/mosaic/gradient_c2l.c
    ${FMS_SRC}/mosaic/interp.c
    ${FMS_SRC}/mosaic/mosaic_util.c
    ${FMS_SRC}/mosaic/read_mosaic.c
    ${FMS_SRC}/affinity/affinity.c
    ${FMS_SRC}/mpp/nsclock.c
    ${FMS_SRC}/mpp/threadloc.c
    ${FMS_SRC}/mpp/mpp_memuse.c
)

set_target_properties(fms PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${CURR_DIR}/lib)
set_target_properties(fms PROPERTIES Fortran_MODULE_DIRECTORY ${CURR_DIR}/mod)

target_include_directories(fms PUBLIC ${FMS_SRC}/include ${FMS_SRC}/fms ${FMS_SRC}/mpp/include ${FMS_SRC}/fms2_io/include )
target_include_directories(fms PRIVATE ${NETCDF_INC_DIR})


