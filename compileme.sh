#!/bin/bash

export FC=mpiifort
export CC=mpiicc
export NetCDF_ROOT=`nc-config --prefix`
export FFLAGS="-O0"

mkdir -p build && cd build
cmake -DCMAKE_INSTALL_PREFIX=$PWD/test_install .. #-DDEBUG=ON ..
make -j4 VERBOSE=1
make install
