#!/bin/bash

FMS_dir="/lustre/f2/dev/Mikyung.Lee/am5_phys/FMS/"
netcdf="`nc-config --fflags` `nc-config --cflags`"

export FC=ftn
export CC=cc

#FMS="-I$FMS_dir/.mods -L$FMS_dir/.libs -I$FMS_dir/include -L$FMS_dir"

export CPPFLAGS="-DINTERNAL_FILE_NML -DCLUBB -Duse_netCDF -DHAVE_SCHED_GETAFFINITY -Duse_LARGEFILE -D__IFC"
export FCFLAGS="-O0 -g -debug -Fno-alias -auto -safe-cray-ptr -ftz -assume byterecl -i4 -r8 -nowarn -sox -traceback -qopenmp `nf-config --fflags`"
export CFLAGS="-g -O0 -debug `nc-config --cflags`"
#export LDFLAGS=$FMS

