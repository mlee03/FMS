#!/bin/bash

export UFS_DIR="/scratch2/GFDL/gfdlscr/Mikyung.Lee/cmake_fms/UFS-weather-model"
export CURR_DIR=$PWD
export FMS_DIR=$UFS_DIR"/FMS"

#: set up modules and module paths UFS uses
source $UFS_DIR/NEMS/src/conf/module-setup.sh.inc
echo $MODULEPATH
#echo $NCEPLIBS

#module load cmake/3.16.1
#module load intel/19.0.5.281
#module load mvapich2/2.3
module load gnu/9.2.0
module load openmpi/3.1.4
module load netcdf/4.7.2
module load hdf5/1.10.5
module list
export NETCDF="/apps/netcdf/4.7.2/gnu/gcc-9.2.0"

#: copy cmake files over over 
cmake_dir="cmake"
mkdir -p $cmake_dir
cp -rv $UFS_DIR/cmake/* $cmake_dir

#: copy source over
#mkdir -p $CURR_DIR/src
#cp -rv $FMS_DIR/* $CURR_DIR/src/.
git clone --recursive --branch 2020.02-beta1 git@github.com:NOAA-GFDL/FMS.git
mv FMS src

export COMPILER=GNU
export FC=gfortran
export CC=gcc
export CXX=gcc
export CMAKE_Fortran_COMPILER=mpif90
export CMAKE_C_COMPILER=mpicc
export CMAKE_CXX_COMPILER=mpicxx
export CMAKE_Platform=hera.gnu
#export CMAKE_C_COMPILER=${CMAKE_C_COMPILER:-mpicc}
#export CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER:-mpicxx}
#export CMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER:-mpif90}

cmake .
