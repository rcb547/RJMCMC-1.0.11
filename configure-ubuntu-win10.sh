#!/bin/bash

### >> mpicc --showme
###    gcc -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi -I/usr/lib/x86_64-linux-gnu/openmpi/include -pthread -L/usr/lib/x86_64-linux-gnu/openmpi/lib -lmpi 

export LIBRJMCMC_ROOT=$PWD/install-ubuntu-win10
./configure  --prefix=$LIBRJMCMC_ROOT --with-openmpi --with-openmpi-include-path=/usr/lib/x86_64-linux-gnu/openmpi/include --with-openmpi-lib-path=/usr/lib/x86_64-linux-gnu/openmpi/lib
make clean
make
make install
