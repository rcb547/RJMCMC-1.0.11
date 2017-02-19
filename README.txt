
RJMCMC Library

This library contains low-level routines for Reversible Jump
Monte-Carlo Markov Chain regression and forward model problems.

Compiling and Installing
------------------------

The library uses the standard autotools build system, to compile:

> ./configure
> make

And install:

> sudo make install

The library can also be linked with an MPI library for parallel
processing. This has been tested withe OpenMPI on Linux and 
MacOSX but other MPI implementations should/may work.

Template Code
-------------

There is a simple template code that can be copied and modified 
to suit your needs in template/part1d. The code is commented to
show how to use the basic features of the library and serves as
a good starting point.

Documentation
-------------

There is documentation for the routines within the code for the
main routines, see include/rjmcmc/regression.h and 
include/rjmcmc/forwardmodel.h. This is written in doxygen format
and can be generated using the dox file in the doc directory, ie

> cd doc
> doxygen rjmcmc.dox

Then the documentation can be viewed in a web browser by opening
the doc/html/index.html file.

Examples
--------

There are several examples in under the examples directory which 
show how to use different aspects of the code. The examples are
arranged in a directory hierarchy corresponding to the broad
categories of the problems they solve, e.g. examples that solve
1D partitioned regression problems are in the 
"examples/1d/partitioned/regression" directory.

Examples that have an "f" as the last character of the name are
fortran examples, e.g. examples/1d/partitioned/fm/functionfitf.

The examples can be built and run inplace, e.g.

> cd examples/1d/partitioned/regression/multistep
> make
> ./multistep

This will run the regression problem then print out the count of
proposals and acceptances and write several files. These results can
be plotted with Gnuplot or R:

> gnuplot
>> plot "./multistep.mean" with lines, "data.txt"

or

> R
>> d <- read.table("data.txt")
>> mean <- read.table("multistep.mean")
>> plot(d$V1, d$V2)
>> lines(mean)

Python Interface
----------------

There is a python interface to the 1D regression codes available under
the python/swig directory with a tutorial introduction under:

  python/tutorial/single/doc
  python/tutorial/multi/doc

The python interface requires swig to build and the tutorials rely on
latex to create the documents and matplotlib for graphing the results.

For Linux, standard installation should follow this procedure:
> ./configure
> make 
> sudo make install
> cd python/swig
> python setup.py build
> sudo python setup.py install

If you want to install the library into a non standard location, 
then you will need to tell pkg-config where to look for its information,
eg for installing to /opt/rjmcmc:

> ./configure --prefix=/opt/rjmcmc
> make
> sudo make install
> export PKG_CONFIG_PATH=/opt/rjmcmc/lib/pkg-config
> cd python/swig
> python setup.py build
> sudo python setup.py install

For MacOSX, the python interface uses the script in macosx/buidlforpython.sh.
So the procedure is:

> make distclean 
> cd macosx
> sh buildforpython.sh
> cd ../python/swig
> python setup.py build
> sudo python setup.py install

Note that the first "make distclean" is only needed if you have previously
configured the source. It will fail harmlessly if you haven't.

Tests
-----

There are some unit tests in the tests directory. These require CUnit to 
be installed to compile.

There are also a series of test programs under the tests/regression
directory that are used to ensure the different variants of functions
are working correctly.

NCI Raijin Compilation
----------------------

On the NCI Raijin supercomputer, using the default OpenMPI module is 
recommended. The following sequence of commands will install the RJMCMC
library to the $HOME/install directory. This can be modified as required.

> cd <path to RJMCMC>/RJMCMC-xx.yy.zz
> module load openmpi
> ./configure --prefix=$HOME/install
> make
> make install

Applications that use the RJMCMC library can then be installed as follows:

> cd <path to APP>/APP-xx.yy.zz
> export PKG_CONFIG_PATH=$HOME/install/lib/pkgconfig
> ./configure --prefix=$HOME/install
> make
> make install

Terrawulf Compilation
---------------------

As of December 2013, the available compilers on the Terrawulf are too
old to compile the Fortran programs that rely on this library as these
rely on some more advanced features of the Fortran standard. Therefore
it is recommended to install locally a more recent GNU Fortran
compiler as follows:

- Download http://gfortran.com/download/x86_64/gcc-4.8-infrastructure.tar.xz
- Download http://gfortran.com/download/x86_64/snapshots/gcc-4.8.tar.xz

Here we assume that all binaries/libraries will be placed under the $HOME/install 
directory. You can modify this as required.

> mkdir -p $HOME/install
> cd $HOME/install
> tar --xz -xf <path to downloads>/gcc-4.8-infrastructure.tar.xz
> tar --xz -x --strip-components=1 -f <path to downloads>/gcc-4.8.tar.xz

The following lines can either be entered directly or added to your $HOME/.bashrc
file:

> export PATH=$PATH:$HOME/install/bin
> export LD_LIBRARY_PATH=$PATH:$HOME/install/lib:$HOME/install/lib64

Next compile the RJMCMC library using the newly install compiler:

> cd <path to RJMCMC>/RJMCMC-xx.yy.zz
> export CC=$HOME/install/bin/gcc
> module load openmpi_144
> export CPPFLAGS=-I$MPI_DIR/include
> export LDFLAGS=-L$MPI_DIR/lib
> ./configure --prefix=$HOME/install
> make
> make install

Note that for MPI programs you will also need to add the Intel library path 
to you LD_LIBRARY_PATH variable as follows:

> export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/intel/Compiler/11.1/038/lib/intel64/

Advanced/Optional
-----------------

To create a Universal Binary on MacOSX, there are some scripts
in the macosx directory. Depending on your version of MacOSX
you may need to edit this, but the process should be:

> make distclean
> cd macosx
> sh builduniversal.sh
> sudo sh installuniversal.sh

Note that the "make distclean" is only required if you have 
already run ./configure after extracting the source. The install
script will install everything under /usr/local. Edit the 
script to modify the install location.

