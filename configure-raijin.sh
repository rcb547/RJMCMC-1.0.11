#!/bin/bash

#On raijin.nci.org.au I changed line 1925 of /short/public/rcb547/apps/iearth/RJMCMC-1.0.11/configure 
#from <am__api_version='1.14'> to <am__api_version='1.11'>
#because on there is a /usr/bin/aclocal-1.11 but not a /usr/bin/aclocal-1.14
	

#module load librjmcmc/gnu
module load librjmcmc/intel
./configure  --with-openmpi --prefix=$LIBRJMCMC_ROOT
make
make install

