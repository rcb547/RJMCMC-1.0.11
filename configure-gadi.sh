#!/bin/bash

#Compilation script for gadi.nci.org.au
module load librjmcmc/gnu
./configure  --with-openmpi --prefix=$LIBRJMCMC_ROOT
make clean
make
make install

#module load librjmcmc/intel
#./configure  CC=icc --with-openmpi --prefix=$LIBRJMCMC_ROOT
#make clean
#make
#make install

