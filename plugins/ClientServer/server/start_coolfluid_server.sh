#!/bin/bash

# export PATH=/software/ar/coolfluid/profiles/fedora9/x86_64/gcc-4.3.0/bin:/software/ar/coolfluid/profiles/fedora9/x86_64/gcc-4.3.0/mpi/openmpi-1.3/bin:$HOME/bin:$PATH
export LD_LIBRARY_PATH=@COOLFluiD_DSO_DIR@:$LD_LIBRARY_PATH

cd @CMAKE_CURRENT_BINARY_DIR@
./app_server $*
