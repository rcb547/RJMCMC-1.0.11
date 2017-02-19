
#
# If you installed rjmcmc to a different location you may have
# to modify these paths
#
INCLUDES = -I/usr/local/include
LIBPATH = -L/usr/lib/openmpi -L/usr/local/lib

#
# If mpi was enabled, you may need an extra include path for mpi.h
#
INCLUDES += -I/usr/include/openmpi

FCC=gfortran
FCFLAGS=-c -g -Wall $(INCLUDES)

LIBS = -lrjmcmc

#
# If mpi was enabled, you will need to include the mpi libraries
#
LIBS += -lmpi

TARGETS = part1d_hierarchical_fortran \
	svdtest

LAPACKLIBS = lapack/liblapack.a blas/libblas.a 

all : $(TARGETS)

lapack/liblapack.a :
	make -C lapack

blas/libblas.a :
	make -C blas

part1d_hierarchical_fortran : part1d_hierarchical_fortran.o svdbacksub.o $(LAPACKLIBS)
	$(FCC) -o part1d_hierarchical_fortran part1d_hierarchical_fortran.o svdbacksub.o $(LIBPATH) $(LAPACKLIBS) $(LIBS)

svdtest : svdbacksub.o svdtest.o $(LAPACKLIBS)
	$(FCC) -o svdtest svdtest.o svdbacksub.o $(LIBPATH) $(LAPACKLIBS) $(LIBS)

%.o : %.f90
	$(FCC) $(FCFLAGS) -o $*.o $*.f90

clean :
	make -C lapack clean
	make -C blas clean
	rm -f $(TARGETS) *.o





