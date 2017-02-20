#!/bin/bash

#On raijin.nci.org.au I changed line 1925 of /short/public/rcb547/apps/iearth/RJMCMC-1.0.11/configure 
#from <am__api_version='1.14'> to <am__api_version='1.11'>
#because on there is a /usr/bin/aclocal-1.11 but not a /usr/bin/aclocal-1.14
	
module load gcc/5.2.0
module load openmpi

./configure --prefix=$PWD/installed/gcc-5.2.0 --with-openmpi
make
make install

